import pandas as pd

input_file = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/reports/ichorCNA_summary.tsv"
output_file = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/reports/ichorCNA_group_summary.tsv"

# Define group membership based on sample ID
sample_to_group = {
    "SHi26-1":  "group1",
    "SHi26-2":  "group1",
    "SHi26-3":  "group2",
    "SHi26-4":  "group2",
    "SHi26-5":  "group3",
    "SHi26-6":  "group4",
    "SHi26-7":  "group5",
    "SHi26-8":  "group5",
    "SHi26-9":  "group6",
    "SHi26-10": "group6",
    "SHi26-11": "group6",
    "SHi26-12": "group6",
}

# Load the summary file
df = pd.read_csv(input_file, sep="\t")

# Add group info
df["Group"] = df["Sample"].map(sample_to_group)

# Group and average
group_summary = df.groupby("Group")[["TumorFraction", "Ploidy"]].mean().round(4)

# Write output
group_summary.to_csv(output_file, sep="\t")

print(f"✅ Group summary written to: {output_file}")

