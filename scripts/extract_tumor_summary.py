import os
import glob
import re

input_dir = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ichor_out"
output_file = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/reports/ichorCNA_summary.tsv"

param_files = glob.glob(os.path.join(input_dir, "SHi26-*.params.txt"))
summary_rows = []

for file in param_files:
    sample = os.path.basename(file).split(".")[0]
    with open(file) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("Tumor Fraction:"):
                tf = line.strip().split(":")[1].strip()
            elif line.startswith("Ploidy:"):
                ploidy = line.strip().split(":")[1].strip()
        summary_rows.append((sample, tf, ploidy))

# Sort numerically by the integer part of SHi26-XX
def sample_sort_key(row):
    match = re.search(r"SHi26-(\d+)", row[0])
    return int(match.group(1)) if match else 0

summary_rows.sort(key=sample_sort_key)

# Write output
with open(output_file, "w") as out:
    out.write("Sample\tTumorFraction\tPloidy\n")
    for row in summary_rows:
        out.write("\t".join(row) + "\n")

print(f"✅ Sorted summary written to: {output_file}")

