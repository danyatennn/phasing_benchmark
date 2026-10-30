#!/usr/bin/env python3
"""Fix two bugs in MethPhaser's meth_phaser_post_processing (get_altered_vcf),
both where the PS of a record in a merged block is rewritten:

1. get_altered_block_start_loc returns the int -1 when the record's phase set
   lies in no merged block, and str.replace() fails on it (TypeError). The
   phase set is now kept as it is.
2. The new PS went in by replace() on the whole sample column. When PS is the
   last FORMAT field (GT:PS, as whatshap writes it) the old value carries the
   line break, replace() took it away, and the record ran into the next one.
   Now only the PS field is set, and the line break stays.

Nothing else changes: which blocks are merged and which are flipped is
MethPhaser's own. The original is kept as <file>.orig; a patched file is left
alone.

  python3 bin/patch_methphaser.py software/methphaser/meth_phaser_post_processing
"""
import re
import shutil
import sys

HELPER = '''

def set_phase_set(sample_column, ps_index, final_blocks):
    # patched (bin/patch_methphaser.py): only the PS field of the sample column
    # is set, the line break it may end with stays, and a phase set that lies
    # in no merged block (get_altered_block_start_loc gives -1) is kept
    body = sample_column.rstrip("\\n")
    fields = body.split(":")
    if ps_index < len(fields) and fields[ps_index].isdigit():
        new = get_altered_block_start_loc(final_blocks, int(fields[ps_index]))
        if new != -1:
            fields[ps_index] = new
    return ":".join(fields) + sample_column[len(body):]
'''

# the function the helper goes after, up to its last line
FUNCTION = re.compile(
    r"^def get_altered_block_start_loc\(final_blocks, original_block_start_loc\):\n"
    r"(?:[ \t]+.*\n)*?[ \t]+return -1\n",
    re.M)

# both branches ("1|0" and "0|1") of get_altered_vcf
CALL = re.compile(
    r"split_rec\[-1\] = split_rec\[-1\]\.replace\(start_loc,\s*"
    r"get_altered_block_start_loc\(current_chrom_final_block,\s*int\(start_loc\)\)\)")
NEW_CALL = "split_rec[-1] = set_phase_set(split_rec[-1], ps_tag_location, current_chrom_final_block)"


def main(path):
    with open(path) as f:
        text = f.read()
    if re.search(r"^def set_phase_set\(", text, re.M):
        print(f"{path}: already patched")
        return 0

    functions = FUNCTION.findall(text)
    calls = CALL.findall(text)
    if len(functions) != 1 or len(calls) != 2:
        print(f"{path}: expected 1 get_altered_block_start_loc and 2 PS rewrites, "
              f"found {len(functions)} and {len(calls)}; not changed", file=sys.stderr)
        return 1

    patched = FUNCTION.sub(lambda m: m.group(0) + HELPER, text, count=1)
    patched = CALL.sub(NEW_CALL, patched)
    compile(patched, path, "exec")

    shutil.copy2(path, path + ".orig")
    with open(path, "w") as f:
        f.write(patched)
    print(f"{path}: patched, the original is {path}.orig")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))
