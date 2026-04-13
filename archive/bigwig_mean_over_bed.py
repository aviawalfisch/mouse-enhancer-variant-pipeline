import sys
import pandas as pd
import pyBigWig

if len(sys.argv) != 4:
    print("Usage: python bigwig_mean_over_bed.py <input.bw> <input.bed> <output.tsv>")
    sys.exit(1)

bw_path = sys.argv[1]
bed_path = sys.argv[2]
out_path = sys.argv[3]

bw = pyBigWig.open(bw_path)

rows = []
with open(bed_path) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        name = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"

        try:
            mean_val = bw.stats(chrom, start, end, type="mean")[0]
            if mean_val is None:
                mean_val = 0.0
        except RuntimeError:
            mean_val = 0.0

        rows.append([name, chrom, start, end, mean_val])

bw.close()

df = pd.DataFrame(rows, columns=["enhancer_id", "chr", "start", "end", "mean_signal"])
df.to_csv(out_path, sep="\t", index=False)

print(f"Saved {len(df)} rows to {out_path}")