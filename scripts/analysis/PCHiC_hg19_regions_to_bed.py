import pandas as pd
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Convert PCHiC regions to bed files')
parser.add_argument('--PC_HiC', type=str, help='Path to PCHiC file')
parser.add_argument('--bed_enhancers', type=str, help='Path to save enhancers bed file')
parser.add_argument('--bed_promoters', type=str, help='Path to save promoters bed file')

args = parser.parse_args()

# Import data
pchic = pd.read_csv(args.PC_HiC, sep="\t")

# Format names
unique_promoters = 'chr'+pchic["baitChr"].astype(str)+'-'+pchic["baitStart"].astype(str)+'-'+pchic["baitEnd"].astype(str)
unique_promoters = pd.Series(unique_promoters.unique())

unique_enhancers = 'chr'+pchic["oeChr"].astype(str)+'-'+pchic["oeStart"].astype(str)+'-'+pchic["oeEnd"].astype(str)
unique_enhancers = pd.Series(unique_enhancers.unique())

# NOT keep old index
unique_enhancers = unique_enhancers.str.split('-', expand=True).reset_index(drop=True)
unique_promoters = unique_promoters.str.split('-', expand=True).reset_index(drop=True)
# Add unique IDs to the peaks
unique_enhancers["ID"] = "enhancer_"+unique_enhancers.index.astype(str)
unique_promoters["ID"] = "promoter_"+unique_promoters.index.astype(str)

# Save bed files
unique_enhancers.to_csv(args.bed_enhancers, sep='\t', index=False, header=False)
unique_promoters.to_csv(args.bed_promoters, sep='\t', index=False, header=False)
