#!/usr/bin/env bash
# Consensus parental-haplotype truth for the reads of a cell: a read is PAT or
# MAT when at least two of three independent-ish labels say so. "Q100" below is
# the sample's diploid assembly (Q100_FASTA, config.sh): Q100 for HG002, the
# HPRC trio assemblies for HG005r2.
#
#   matpat   the Q100 copy the read aligns to better, alignment-score margin
#            (bin/matpat_truth.sh, reads/<cell>.matpat.tsv.gz)
#   kmer     the share of the read's parent-specific 31-mers (bin/yak_truth.sh,
#            reads/<cell>.kmer.tsv.gz)
#   q100snp  the alleles the read carries at the SNPs where the two Q100 copies
#            differ (this script, reads/<cell>.q100snp.tsv.gz)
#
# q100snp, made here in three steps:
#   1. once: the paternal copy aligned to the maternal one (minimap2 -cx asm5
#      --cs), paftools.js call, and the informative SNPs kept: one alignment
#      over the site, mapQ >= 30, biallelic A/C/G/T, the two alleles different
#        <Q100_WORK>/informative_snps.tsv   contig 0-based-pos maternal paternal
#   2. per cell: the reads (primary records of the cell's BAM) aligned to the
#      maternal copy (lr:hq for ONT, map-hifi for PacBio) and streamed, never
#      stored. A read counts its maternal and paternal alleles at those SNPs,
#      base quality >= Q100SNP_MIN_BASEQ; reads under Q100SNP_MIN_MAPQ, and
#      with Q100SNP_EXCLUDE_SPLIT=1 reads with a supplementary alignment, are
#      not labelled. Of the loci where the read shows one of the two alleles,
#      at least Q100SNP_MIN_LOCI, a share of at least Q100SNP_MIN_RATIO for one
#      parent: p / m; enough loci but no clear parent: a; too few: 0
#   3. per cell: the consensus of matpat, kmer and q100snp
#        reads/<cell>.consensus.tsv.gz   read label matpat kmer q100snp agree
#      CONSENSUS_RULE  2of3 (default)   PAT if >= 2 say PAT, the third may disagree
#                      2of3_strict      >= 2 agree and none says the other parent
#                      all3             all three agree
#      The rule is written in the header; a table made with another rule is
#      rebuilt, the q100snp tables are kept.
#
# Needs per cell: its BAM, reads/<cell>.matpat.tsv.gz and reads/<cell>.kmer.tsv.gz.
# Then READ_TRUTH=consensus bash bin/read_recall_strata.sh, as for matpat.
#
#   PHASING_DEPTHS=30 REFERENCE_NAME=T2T bash bin/consensus_truth.sh
#   CONSENSUS_RULE=2of3_strict PHASING_DEPTHS=30 bash bin/consensus_truth.sh

set -euo pipefail

source "/autofs/bal19/zxzheng/somatic/Clair-somatic/scripts/data_config_germline.sh"
source "/autofs/bal36/dten/phasing_benchmark/lib.sh"
source "/autofs/bal36/dten/phasing_benchmark/config.sh"

[[ -n "${Q100_FASTA}" ]] || die "no diploid assembly for ${SAMPLE}: set diploid_fasta in ${SAMPLES_TSV} (bin/asm_truth.sh writes one)"
require_file "${Q100_FASTA}" "diploid HG002 Q100 assembly"

CONSENSUS_THREADS="${CONSENSUS_THREADS:-32}"
CONSENSUS_RULE="${CONSENSUS_RULE:-2of3}"
CONSENSUS_PYTHON="${CONSENSUS_PYTHON:-/autofs/bal36/dten/miniforge3/envs/phasing/bin/python}"   # has pysam
Q100_WORK="${Q100_WORK:-$(dirname "${Q100_FASTA}")/informative}"
Q100SNP_MIN_LOCI="${Q100SNP_MIN_LOCI:-5}"
Q100SNP_MIN_RATIO="${Q100SNP_MIN_RATIO:-0.90}"
Q100SNP_MIN_MAPQ="${Q100SNP_MIN_MAPQ:-20}"
Q100SNP_MIN_BASEQ="${Q100SNP_MIN_BASEQ:-20}"
Q100SNP_EXCLUDE_SPLIT="${Q100SNP_EXCLUDE_SPLIT:-1}"
# paftools.js and its k8 come with the minimap2 binary release
K8="${K8:-$(dirname "${MINIMAP2}")/k8}"
PAFTOOLS="${PAFTOOLS:-$(dirname "${MINIMAP2}")/paftools.js}"

case "${CONSENSUS_RULE}" in
  2of3|2of3_strict|all3) ;;
  *) die "CONSENSUS_RULE must be 2of3, 2of3_strict or all3, not ${CONSENSUS_RULE}" ;;
esac
require_file "${K8}" "k8 of the minimap2 release (set K8)"
require_file "${PAFTOOLS}" "paftools.js of the minimap2 release (set PAFTOOLS)"
"${CONSENSUS_PYTHON}" -c 'import pysam' 2>/dev/null || die "${CONSENSUS_PYTHON} has no pysam (set CONSENSUS_PYTHON)"
mkdir -p "${Q100_WORK}" "${READS_DIR}" "${LOGS_DIR}/reads" "${TMP_DIR}"

SNPS="${Q100_WORK}/informative_snps.tsv"

# the two steps that read one record at a time: assign (SAM on stdin -> one
# line per read) and merge (three label tables -> the consensus)
IFS= read -r -d '' PYCODE <<'PY' || true
import bisect, gzip, sys

def opener(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def assign(args):
    import pysam
    snps, out = args[0], args[1]
    min_loci, min_ratio, min_mapq, min_baseq, exclude_split = (
        int(args[2]), float(args[3]), int(args[4]), int(args[5]), args[6] == "1")
    pos, mat, pat = {}, {}, {}
    with open(snps) as fh:
        for line in fh:
            c, p, m, q = line.rstrip("\n").split("\t")
            pos.setdefault(c, []).append(int(p))
            mat.setdefault(c, []).append(m)
            pat.setdefault(c, []).append(q)
    n = {"p": 0, "m": 0, "a": 0, "0": 0}
    sam = pysam.AlignmentFile("-", "r")
    with gzip.open(out, "wt", compresslevel=3) as fo:
        fo.write("#read\tlabel\tclass\tpat\tmat\tother\tmapq\tsplit\n")
        for r in sam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            split = int(r.has_tag("SA"))
            p_ = m_ = o_ = 0
            c = r.reference_name
            if c in pos and r.mapping_quality >= min_mapq and not (exclude_split and split):
                P, M, Q = pos[c], mat[c], pat[c]
                k = bisect.bisect_left(P, r.reference_start)
                j = bisect.bisect_left(P, r.reference_end)
                if k < j:
                    seq, qual = r.query_sequence, r.query_qualities
                    rp, qp = r.reference_start, 0
                    for op, ln in r.cigartuples:
                        if op in (0, 7, 8):                     # M = X
                            while k < j and P[k] < rp + ln:
                                i = qp + P[k] - rp
                                if qual is None or qual[i] >= min_baseq:
                                    b = seq[i].upper()
                                    if b == M[k]:
                                        m_ += 1
                                    elif b == Q[k]:
                                        p_ += 1
                                    else:
                                        o_ += 1
                                k += 1
                            rp += ln
                            qp += ln
                        elif op in (1, 4):                      # I S
                            qp += ln
                        elif op in (2, 3):                      # D N: no base there
                            while k < j and P[k] < rp + ln:
                                k += 1
                            rp += ln
                        if k >= j:
                            break
            t = p_ + m_
            if t < min_loci:
                cls = "0"
            elif p_ >= min_ratio * t:
                cls = "p"
            elif m_ >= min_ratio * t:
                cls = "m"
            else:
                cls = "a"
            n[cls] += 1
            label = {"p": "PAT", "m": "MAT"}.get(cls, "none")
            fo.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n" % (r.query_name, label, cls, p_, m_, o_,
                                                          r.mapping_quality, split))
    print("%s: p %d, m %d, a %d, 0 %d" % (out, n["p"], n["m"], n["a"], n["0"]), file=sys.stderr)

def labels(path):
    out = {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 1 and f[1] in ("PAT", "MAT") and f[0] not in out:
                out[f[0]] = f[1]
    return out

def merge(args):
    rule, out, paths = args[0], args[1], args[2:5]
    names = ["matpat", "kmer", "q100snp"]
    tabs = [labels(p) for p in paths]
    reads = set().union(*tabs)
    pattern = {}
    n = {"PAT": 0, "MAT": 0, "none": 0}
    with gzip.open(out, "wt", compresslevel=3) as fo:
        fo.write("#read\tlabel\t%s\tagree\trule=%s\n" % ("\t".join(names), rule))
        for r in sorted(reads):
            v = [t.get(r, ".") for t in tabs]
            np_, nm = v.count("PAT"), v.count("MAT")
            if rule == "2of3":
                label = "PAT" if np_ >= 2 else "MAT" if nm >= 2 else "none"
            elif rule == "2of3_strict":
                label = "PAT" if np_ >= 2 and nm == 0 else "MAT" if nm >= 2 and np_ == 0 else "none"
            else:
                label = "PAT" if np_ == 3 else "MAT" if nm == 3 else "none"
            agree = "%d/%d" % (max(np_, nm), min(np_, nm))
            n[label] += 1
            pattern[agree] = pattern.get(agree, 0) + 1
            fo.write("%s\t%s\t%s\t%s\n" % (r, label, "\t".join(v), agree))
    print("%s: rule %s, PAT %d, MAT %d, none %d; labels agreeing/opposing: %s" % (
        out, rule, n["PAT"], n["MAT"], n["none"],
        ", ".join("%s %d" % (k, pattern[k]) for k in sorted(pattern, reverse=True))), file=sys.stderr)
    for name, t in zip(names, tabs):
        print("  %s labels %d" % (name, len(t)), file=sys.stderr)

{"assign": assign, "merge": merge}[sys.argv[1]](sys.argv[2:])
PY

# ----------------------------------------------- 1. informative Q100 SNPs

if newer_than "${SNPS}" "${Q100_FASTA}"; then
  log "informative Q100 SNPs: $(wc -l <"${SNPS}") (kept)"
else
  logf="${LOGS_DIR}/reads/q100_informative.log"
  if [[ ! -s "${Q100_WORK}/maternal.fa" || ! -s "${Q100_WORK}/paternal.fa" ]]; then
    log "splitting $(basename "${Q100_FASTA}") into its maternal and paternal copies"
    gzip -cdf "${Q100_FASTA}" | awk -v mat="${Q100_WORK}/maternal.fa.tmp" -v pat="${Q100_WORK}/paternal.fa.tmp" '
      /^>/ { h = toupper(substr($1, 2))
             w = (h ~ /[_.#-](MATERNAL|MAT)$/) ? 1 : (h ~ /[_.#-](PATERNAL|PAT)$/) ? 2 : 0 }
      w == 1 { print > mat }
      w == 2 { print > pat }'
    [[ -s "${Q100_WORK}/maternal.fa.tmp" && -s "${Q100_WORK}/paternal.fa.tmp" ]] \
      || die "no contig named as maternal or paternal in ${Q100_FASTA}"
    mv -f "${Q100_WORK}/maternal.fa.tmp" "${Q100_WORK}/maternal.fa"
    mv -f "${Q100_WORK}/paternal.fa.tmp" "${Q100_WORK}/paternal.fa"
  fi
  log "paternal copy -> maternal copy (minimap2 asm5)"
  "${MINIMAP2}" -cx asm5 --cs -t "${CONSENSUS_THREADS}" \
      "${Q100_WORK}/maternal.fa" "${Q100_WORK}/paternal.fa" 2>"${logf}" \
    | sort -k6,6 -k8,8n -S 8G -T "${TMP_DIR}" >"${Q100_WORK}/pat_to_mat.sorted.paf.tmp" \
    || die "pat-to-mat alignment failed, see ${logf}"
  mv -f "${Q100_WORK}/pat_to_mat.sorted.paf.tmp" "${Q100_WORK}/pat_to_mat.sorted.paf"
  log "paftools.js call"
  "${K8}" "${PAFTOOLS}" call "${Q100_WORK}/pat_to_mat.sorted.paf" \
      >"${Q100_WORK}/pat_to_mat.variants.txt.tmp" 2>>"${logf}" \
    || die "paftools.js call failed, see ${logf}"
  mv -f "${Q100_WORK}/pat_to_mat.variants.txt.tmp" "${Q100_WORK}/pat_to_mat.variants.txt"
  # V ctg start end depth mapq ref(maternal) alt(paternal) ...
  awk 'BEGIN { OFS = "\t" }
    $1 == "V" && $5 == 1 && $6 >= 30 && length($7) == 1 && length($8) == 1 &&
    toupper($7) ~ /^[ACGT]$/ && toupper($8) ~ /^[ACGT]$/ && toupper($7) != toupper($8) {
      print $2, $3, toupper($7), toupper($8) }' "${Q100_WORK}/pat_to_mat.variants.txt" \
    | sort -k1,1 -k2,2n -u -T "${TMP_DIR}" >"${SNPS}.tmp"
  [[ -s "${SNPS}.tmp" ]] || die "no informative SNP kept, see ${Q100_WORK}/pat_to_mat.variants.txt"
  mv -f "${SNPS}.tmp" "${SNPS}"
  log "informative Q100 SNPs: $(wc -l <"${SNPS}") -> ${SNPS}"
fi

preset_of() {
  case "$1" in
    ont4k|ont5k) printf 'lr:hq\n' ;;
    hifi|revio)  printf 'map-hifi\n' ;;
  esac
}

# ------------------------------------------------- 2. and 3. per cell

for platform in "${PLATFORMS[@]}"; do
  for depth in "${DEPTHS[@]}"; do
    cell="$(bam_key "${platform}" "${REFERENCE_NAME}" "${depth}")"
    bam="$(sub_bam "${platform}" "${REFERENCE_NAME}" "${depth}")"
    [[ -s "${bam}" ]] || { warn "no BAM for ${cell}"; continue; }
    snp_out="${READS_DIR}/${cell}.q100snp.tsv.gz"
    logf="${LOGS_DIR}/reads/${cell}.q100snp.log"

    if newer_than "${snp_out}" "${bam}" && newer_than "${snp_out}" "${SNPS}"; then
      log "${cell}: q100snp kept"
    else
      preset="$(preset_of "${platform}")"
      mmi="${Q100_WORK}/maternal.${preset//:/_}.mmi"
      if [[ ! -s "${mmi}" ]]; then
        log "minimap2 index of the maternal copy for ${preset}"
        "${MINIMAP2}" -t "${CONSENSUS_THREADS}" -x "${preset}" -d "${mmi}.tmp" "${Q100_WORK}/maternal.fa"
        mv -f "${mmi}.tmp" "${mmi}"
      fi
      log "${cell}: reads -> maternal Q100 copy (${preset}), alleles at the informative SNPs"
      "${SAMTOOLS}" fastq -@ 4 -n -F 0x900 "${bam}" 2>"${logf}" \
        | "${MINIMAP2}" -t "${CONSENSUS_THREADS}" -ax "${preset}" --secondary=no "${mmi}" - 2>>"${logf}" \
        | "${CONSENSUS_PYTHON}" -c "${PYCODE}" assign "${SNPS}" "${snp_out}.tmp.gz" \
            "${Q100SNP_MIN_LOCI}" "${Q100SNP_MIN_RATIO}" "${Q100SNP_MIN_MAPQ}" \
            "${Q100SNP_MIN_BASEQ}" "${Q100SNP_EXCLUDE_SPLIT}" 2>>"${logf}" \
        || die "q100snp failed for ${cell}, see ${logf}"
      mv -f "${snp_out}.tmp.gz" "${snp_out}"
      tail -1 "${logf}" >&2
    fi

    matpat="${READS_DIR}/${cell}.matpat.tsv.gz"
    kmer="${READS_DIR}/${cell}.kmer.tsv.gz"
    out="${READS_DIR}/${cell}.consensus.tsv.gz"
    [[ -s "${matpat}" ]] || { warn "${cell}: no ${matpat} (bin/matpat_truth.sh), no consensus"; continue; }
    [[ -s "${kmer}" ]] || { warn "${cell}: no ${kmer} (bin/yak_truth.sh), no consensus"; continue; }
    rule_now="$(gzip -dc "${out}" 2>/dev/null | head -1 || true)"
    if newer_than "${out}" "${matpat}" && newer_than "${out}" "${kmer}" && newer_than "${out}" "${snp_out}" \
        && [[ "${rule_now}" == *"rule=${CONSENSUS_RULE}"* ]]; then
      log "${cell}: consensus kept"
      continue
    fi
    "${CONSENSUS_PYTHON}" -c "${PYCODE}" merge "${CONSENSUS_RULE}" "${out}.tmp.gz" \
        "${matpat}" "${kmer}" "${snp_out}" 2>>"${logf}" \
      || die "consensus failed for ${cell}, see ${logf}"
    mv -f "${out}.tmp.gz" "${out}"
    tail -4 "${logf}" >&2
  done
done
log "consensus truth: ${READS_DIR}/*.consensus.tsv.gz"
