#!/usr/bin/env python3
"""A Word document with the figures, a caption under each and the metric
definitions, for sharing drafts: built with pandoc from the PNGs that
make_plots.py --png wrote, so the document shows exactly what the code draws.

  python plots/make_report.py --figures plots/pdf_no_margin --out report.docx
  python plots/make_report.py ... --phasers "WhatsHap, LongPhase and HapCUT2"

A figure is included when its PNG is there (both references for the
per-reference ones, whichever exist); the rest are left out without a gap in
the numbering. Needs pandoc on the PATH or --pandoc.
"""

import argparse
import datetime
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

REFS = [("GRCh38", "GRCh38"), ("T2T", "T2T-CHM13")]

TRUTH_INPUT = ("Phasers on the GIAB truth input (the unphased truth het SNPs), so only the "
               "phaser differs.")

# (file name prefix, per reference?, title, caption)
FIGURES = [
    ("fig1a_switch_rate_by_platform", True, "Switch error rate against coverage",
     TRUTH_INPUT + " Log scale; lower is better."),
    ("fig1b_block_n50_by_platform", True, "Phase block NG50 against coverage",
     "Same runs as above. NG50 over the length of chr1-22. Log scale; higher is better."),
    ("fig1c_truth_phased_by_platform", True, "Truth het SNPs phased",
     "Same runs. Share of the benchmark het SNPs that are in a phase block of two or more SNPs."),
    ("fig3_accuracy_vs_contiguity", True, "Accuracy against contiguity",
     TRUTH_INPUT + " One line per phaser, from 10x (small marker) to 50x (large marker). The "
     "switch error rate axis is inverted, so the top right corner is best."),
    ("fig4_switch_vs_flip_30x", True, "Switch and flip errors at 30x",
     TRUTH_INPUT + " Dark: switch errors (a stretch phased the wrong way). Light: flip errors "
     "(a single site phased the wrong way). Per 100 assessed pairs of SNPs."),
    ("figS1_switch_rate_per_chromosome_30x", True, "Switch error rate per autosome at 30x",
     TRUTH_INPUT + " Box: the 22 autosomes; points: single autosomes. Autosomes without a "
     "switch error are drawn on the line labelled 0."),
    ("fig5a_caller_fnr", True, "Variant callers: false negative rate",
     "Biallelic heterozygous PASS SNPs of each caller against the truth, inside the benchmark "
     "regions. Log scale; lower is better."),
    ("fig5b_caller_fdr", True, "Variant callers: false discovery rate",
     "Same calls. Log scale; lower is better."),
    ("fig5d_truth_recovered_30x", True, "What calling and what phasing lose, 30x",
     "Both ends are shares of the same benchmark het SNPs. Open circle: called as het SNP. "
     "Filled circle: called and phased correctly. The gap is lost in phasing; the distance of "
     "the open circle from 100% is lost in calling. Caller VCFs were given to the phasers as "
     "they are."),
    ("fig5c_truth_lost_all_callers", True, "Truth het SNPs not called or phased wrongly, every run",
     "Rows: input x phaser; columns: platform x coverage. Percent of the benchmark het SNPs, "
     "log colour scale; lower is better."),
    ("figS2a_phasing_lost_all_callers", True, "Loss in phasing, every run",
     "Truth het SNPs the caller did find but that were left unphased or phased wrongly, "
     "percent of the benchmark het SNPs."),
    ("figS2b_block_ng50_all_callers", True, "Phase block NG50, every run",
     "NG50 over the length of chr1-22, Mb."),
    ("fig5e_switch_rate_vs_truth_phased_30x", True, "Switch error rate against the share of the truth phased, 30x",
     "Colour: input; shape: phaser. The switch error rate of a caller is measured on the sites "
     "that caller found, so it is read together with the x axis."),
    ("fig6a_reference_switch_rate_scatter_by_platform", False, "GRCh38 against T2T-CHM13: switch error rate",
     "Same run on both references, truth input, every coverage. Dashed line: equal on both."),
    ("fig6a_reference_switch_rate_ratio_by_technology", False, "GRCh38 against T2T-CHM13: switch error rate ratio",
     "Each point: one run (same platform, coverage and phaser) on both references, truth input; "
     "a box pools the five coverages. Box: interquartile range; black line: median; number "
     "above: the median change on T2T-CHM13 in percent; grey band: within 10%. Reads were aligned "
     "differently on the two references (provider alignments on GRCh38, minimap2 2.31 on "
     "T2T-CHM13)."),
    ("fig6b_reference_block_ng50_scatter_by_platform", False, "GRCh38 against T2T-CHM13: phase block NG50",
     "As above, for NG50 over the length of chr1-22."),
    ("fig6b_reference_block_ng50_ratio_by_technology", False, "GRCh38 against T2T-CHM13: NG50 ratio",
     "As above, for NG50."),
    ("fig6c_reference_recall_phased", False, "GRCh38 against T2T-CHM13: truth het SNPs phased correctly",
     "Same run on both references, truth input. Dashed line: equal on both."),
    ("fig7a_bases_tagged", True, "Read bases assigned a haplotype",
     TRUTH_INPUT + " Share of read bases on chr1-22 in reads tagged in a phase block of two or "
     "more SNPs."),
    ("figS4_reads_tagged", True, "Reads assigned a haplotype (by read count)",
     "As above, counting reads instead of bases. A read without a het SNP cannot be tagged, and "
     "ONT has many short reads, so by read count ONT looks worse than by bases."),
    ("fig7b_read_precision", True, "Read-level precision",
     TRUTH_INPUT + " Reads on the truth haplotype, of the reads tagged both by the run and by "
     "the truth. Phase blocks oriented by their SNPs."),
    ("fig7c_read_recall", True, "Read-level recall",
     TRUTH_INPUT + " Reads on the truth haplotype, of the reads tagged by the truth."),
    ("fig7d_anchor_vs_variant_orientation", True, "Two ways to orient phase blocks",
     "Every run. Reads on the wrong haplotype with blocks oriented by their SNPs (x) and by an "
     "anchor read (y). Dashed line: equal."),
    ("fig8_switch_rate_by_stratum_30x", True, "Switch error rate by genomic context, 30x",
     TRUTH_INPUT + " GIAB stratifications v3.6, inside the benchmark regions, chr1-22."),
]

DEFINITIONS = r"""
# Metrics

*T*: the benchmark het SNPs: phased, biallelic, heterozygous SNPs of the GIAB truth inside its benchmark regions, chr1-22. Every phasing metric is computed on these sites.

**Switch error rate.** *S*: switch errors; *P*: assessed pairs, i.e. pairs of consecutive benchmark het SNPs phased in the same block by the run (whatshap compare).

$$\mathrm{SER} = \frac{S}{P} \times 100\%$$

**Switch and flip errors.** A flip is a single site phased the wrong way (two switches in a row) and is counted once; the remaining switches are switch errors. Both per 100 assessed pairs.

**Truth het SNPs phased.** *C*: benchmark het SNPs in a phase block of two or more phased SNPs.

$$\mathrm{Phased} = \frac{C}{|T|} \times 100\%$$

**Truth het SNPs phased correctly.** *H*: blockwise Hamming distance, the sites to flip in each block to match the truth.

$$\mathrm{Correct} = \frac{C - H}{|T|} \times 100\%$$

**Phase block NG50.** Blocks sorted by span, *L*~1~ ≥ *L*~2~ ≥ …; *G*: total length of chr1-22.

$$\mathrm{NG50} = L_k, \quad k = \min\Big\{ j : \sum_{i=1}^{j} L_i \ge G/2 \Big\}$$

**Variant calling.** *TP*: called het SNPs matching a benchmark het SNP (chromosome, position, REF, ALT); *FP*: other calls inside the benchmark regions; *FN* = |*T*| − *TP*.

$$\mathrm{FNR} = \frac{FN}{|T|}, \qquad \mathrm{FDR} = \frac{FP}{TP + FP}$$

**Lost in calling and in phasing** (shares of |*T*|).

$$\mathrm{Called} = \frac{TP}{|T|}, \qquad \mathrm{Lost\ in\ calling} = 100\% - \mathrm{Called}, \qquad \mathrm{Lost\ in\ phasing} = \mathrm{Called} - \mathrm{Correct}$$

**Truth het SNPs not called or phased wrongly** (every run; the complement of Correct):

$$\mathrm{Lost} = 100\% - \mathrm{Correct} = \Big(1 - \frac{C - H}{|T|}\Big) \times 100\%$$

**Called truth het SNPs left unphased or phased wrongly** (the loss in phasing alone, as a share of the truth):

$$\mathrm{Lost\ in\ phasing} = \frac{TP - (C - H)}{|T|} \times 100\%$$

**Reads and bases assigned a haplotype.** A read counts as tagged if it has a haplotype tag (whatshap haplotag) in a block of two or more phased SNPs. Primary alignments on chr1-22.

$$\mathrm{Reads\ tagged} = \frac{N_{\mathrm{tagged}}}{N_{\mathrm{reads}}}, \qquad \mathrm{Bases\ tagged} = \frac{\sum_{\mathrm{tagged}} \ell}{\sum_{\mathrm{all}} \ell}$$

with *ℓ* the read length.

**Read-level precision and recall.** Haplotype labels are arbitrary within a block, so each block is first oriented: by its SNPs against the truth genotypes, or by an anchor read (the longest read of the block with a truth label). *R*~correct~: reads on the truth haplotype after orientation; *R*~both~: reads tagged by the run and by the truth; *R*~truth~: reads tagged by the truth.

$$\mathrm{Precision} = \frac{R_{\mathrm{correct}}}{R_{\mathrm{both}}}, \qquad \mathrm{Recall} = \frac{R_{\mathrm{correct}}}{R_{\mathrm{truth}}}$$

**GRCh38 against T2T-CHM13.** For the same platform, coverage, input and phaser:

$$\mathrm{Ratio} = \frac{\mathrm{metric\ on\ T2T\text{-}CHM13}}{\mathrm{metric\ on\ GRCh38}}$$
"""


def setup_text(phasers):
    return """# Setup

- Sample: HG002. Truth: GIAB v5.0q small variants (Q100-based) on GRCh38 and on T2T-CHM13v2.0.
- Data: ONT R10.4.1 4 kHz, ONT R10.4.1 5 kHz, PacBio HiFi Sequel II, PacBio Revio, each downsampled to 10, 20, 30, 40 and 50x (samtools view -s, fixed seed).
- Inputs to the phasers: the unphased GIAB truth (phaser alone), Clair3, DeepVariant, longcallD.
- Phasers shown: %s.
- All metrics on chr1-22.
""" % phasers


def arial_theme(path):
    """pandoc's default theme uses Calibri / Cambria; make every font Arial."""
    tmp = path + ".tmp"
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename in ("word/theme/theme1.xml", "word/styles.xml"):
                text = data.decode("utf-8")
                text = re.sub(r'typeface="(Calibri|Cambria|Calibri Light|Aptos|Aptos Display)"',
                              'typeface="Arial"', text)
                text = re.sub(r'w:(ascii|hAnsi|cs|eastAsia)="(Calibri|Cambria|Aptos)"', r'w:\1="Arial"', text)
                if item.filename == "word/styles.xml":
                    # theme fonts are not honoured by every viewer: name Arial outright
                    text = re.sub(r"<w:rFonts [^>]*/>", '<w:rFonts w:ascii="Arial" w:hAnsi="Arial" '
                                  'w:cs="Arial" w:eastAsia="Arial"/>', text)
                data = text.encode("utf-8")
            zout.writestr(item, data)
    os.replace(tmp, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--figures", required=True, help="folder with the PNGs of make_plots.py --png")
    parser.add_argument("--out", required=True, help="the .docx to write")
    parser.add_argument("--phasers", default="WhatsHap, LongPhase and HapCUT2")
    parser.add_argument("--title", default="Phasing benchmark: figures")
    parser.add_argument("--pandoc", default=shutil.which("pandoc") or "pandoc")
    args = parser.parse_args()

    lines = ["---", 'title: "%s"' % args.title,
             'subtitle: "Draft, %s"' % datetime.date.today().isoformat(), "---", "",
             setup_text(args.phasers), DEFINITIONS, "", "# Figures", ""]
    n = 0
    included = []
    for prefix, per_ref, title, caption in FIGURES:
        if per_ref:
            panels = [(label, os.path.join(args.figures, "%s_%s.png" % (prefix, ref)))
                      for ref, label in REFS]
        else:
            panels = [("", os.path.join(args.figures, prefix + ".png"))]
        panels = [(l, p) for l, p in panels if os.path.exists(p)]
        if not panels:
            continue
        n += 1
        included.append(prefix)
        lines.append("## Figure %d. %s" % (n, title))
        lines.append("")
        for i, (label, path) in enumerate(panels):
            alt = "(%s) %s" % ("ab"[i], label) if per_ref and len(panels) > 1 else label
            lines.append("![%s](%s){width=100%%}" % (alt, os.path.abspath(path)))
            lines.append("")
        if per_ref and len(panels) > 1:
            caption += " HG002, (a) GRCh38, (b) T2T-CHM13."
        elif per_ref:
            caption += " HG002, %s." % panels[0][0]
        else:
            caption += " HG002."
        lines.append(caption)
        lines.append("")
        lines.append("\\newpage")
        lines.append("")
    if not n:
        sys.exit("no figure PNG found in %s - run make_plots.py with --png" % args.figures)

    with tempfile.TemporaryDirectory() as tmp:
        md = os.path.join(tmp, "report.md")
        with open(md, "w") as fh:
            fh.write("\n".join(lines))
        out = os.path.abspath(args.out)
        subprocess.run([args.pandoc, md, "-f", "markdown+tex_math_dollars+implicit_figures",
                        "-o", out + ".tmp.docx"], check=True)
        arial_theme(out + ".tmp.docx")
        os.replace(out + ".tmp.docx", out)
    print("%d figures -> %s" % (n, args.out))
    for prefix in included:
        print("  %s" % prefix)


if __name__ == "__main__":
    main()
