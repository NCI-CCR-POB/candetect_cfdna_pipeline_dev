"""Per-library count files from extract_5p_motifs.sh -> combined motif matrices.

  python aggregate_motif_counts.py <counts_dir> <counts.tsv> <fractions.tsv>
"""
import sys
from itertools import product
from pathlib import Path

import pandas as pd

if len(sys.argv) != 4:
    sys.exit(
        "Usage: aggregate_motif_counts.py "
        "<motif_tsv_dir> <combined_counts.tsv> <combined_fractions.tsv>"
    )

MOTIFS = ["".join(p) for p in product("ACGT", repeat=4)]
src = Path(sys.argv[1])
output_counts = Path(sys.argv[2])
output_fractions = Path(sys.argv[3])
output_counts.parent.mkdir(parents=True, exist_ok=True)
output_fractions.parent.mkdir(parents=True, exist_ok=True)

cols, other = {}, {}
for f in sorted(src.glob("*.tsv")):
    d = pd.read_csv(f, sep="\t", header=None, index_col=0)[1]
    cols[f.stem] = d.reindex(MOTIFS).fillna(0).astype(int)
    other[f.stem] = int(d.get("OTHER", 0))
if not cols:
    sys.exit(f"[FAIL] no *.tsv under {src}")
counts = pd.DataFrame(cols).rename_axis("motif")

totals = counts.sum(axis=0)
if (totals == 0).any():
    bad = ", ".join(totals[totals == 0].index)
    sys.exit(f"[FAIL] zero valid ACGT motif counts for: {bad}")

counts.to_csv(output_counts, sep="\t")
fractions = counts.div(totals, axis=1)
fractions.to_csv(output_fractions, sep="\t")

print(f"{len(cols)} libraries; ends per library median {totals.median()/1e6:.1f}M [{totals.min()/1e6:.1f}, {totals.max()/1e6:.1f}]; "
      f"OTHER max {max(other.values())}", file=sys.stderr)
