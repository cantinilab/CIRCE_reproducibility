import pandas as pd
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Convert PCHiC regions to bed files')
parser.add_argument('--PC_HiC', type=str, help='Path to PCHiC file')
parser.add_argument('--hg19_enhancers', type=str, help='Path to save hg19 enhancers bed file')
parser.add_argument('--hg38_enhancers', type=str, help='Path to save hg38 enhancers bed file')
parser.add_argument('--hg19_promoters', type=str, help='Path to save hg19 promoters bed file')
parser.add_argument('--hg38_promoters', type=str, help='Path to save hg38 promoters bed file')
parser.add_argument('--hg38_PCHiC', type=str, help='Path to save PCHiC file with hg38 coordinates')

args = parser.parse_args()

# Import data
pchic = pd.read_csv(args.PC_HiC, sep="\t")

### Enhancers
hg19_enhancers = pd.read_csv(args.hg19_enhancers, sep="\t", header=None).astype(str)
hg38_enhancers = pd.read_csv(args.hg38_enhancers, sep='\t', header=None).astype(str)
hg19_enhancers.columns = ["hg19_chrom", "hg19_start", "hg19_end", "enhancer_key"]
hg38_enhancers.columns = ["hg38_chrom", "hg38_start", "hg38_end", "enhancer_key"]
hg19_enhancers["hg19_region"] = hg19_enhancers["hg19_chrom"]+":"+hg19_enhancers["hg19_start"]+":"+hg19_enhancers["hg19_end"]
hg38_enhancers["hg38_region"] = hg38_enhancers["hg38_chrom"]+":"+hg38_enhancers["hg38_start"]+":"+hg38_enhancers["hg38_end"]

hg19_to_hg38 = hg38_enhancers.merge(hg19_enhancers, on="enhancer_key", how="inner")

# Update PCHiC coordinates for enhancers
pchic["enhancer_hg19"] = "chr"+pchic["oeChr"].astype(str)+':'+pchic["oeStart"].astype(str)+':'+pchic["oeEnd"].astype(str)
pchic["enhancer_hg38"] = pchic["enhancer_hg19"].map(hg19_to_hg38[["hg19_region", "hg38_region"]].set_index("hg19_region")["hg38_region"].to_dict())
pchic[["oeChr", "oeStart", "oeEnd"]] = pchic["enhancer_hg38"].str.split(':', expand=True)

### Promoters
hg19_promoters = pd.read_csv(args.hg19_promoters, sep="\t", header=None).astype(str)
hg38_promoters = pd.read_csv(args.hg38_promoters, sep='\t', header=None).astype(str)
hg19_promoters.columns = ["hg19_chrom", "hg19_start", "hg19_end", "promoter_key"]
hg38_promoters.columns = ["hg38_chrom", "hg38_start", "hg38_end", "promoter_key"]
hg19_promoters["hg19_region"] = hg19_promoters["hg19_chrom"]+":"+hg19_promoters["hg19_start"]+":"+hg19_promoters["hg19_end"]
hg38_promoters["hg38_region"] = hg38_promoters["hg38_chrom"]+":"+hg38_promoters["hg38_start"]+":"+hg38_promoters["hg38_end"]

hg19_to_hg38_promoters = hg38_promoters.merge(hg19_promoters, on="promoter_key", how="inner")

# Update PCHiC coordinates for promoters
pchic["promoter_hg19"] = "chr"+pchic["baitChr"].astype(str)+':'+pchic["baitStart"].astype(str)+':'+pchic["baitEnd"].astype(str)
pchic["promoter_hg38"] = pchic["promoter_hg19"].map(hg19_to_hg38_promoters[["hg19_region", "hg38_region"]].set_index("hg19_region")["hg38_region"].to_dict())
pchic[["baitChr", "baitStart", "baitEnd"]] = pchic["promoter_hg38"].str.split(':', expand=True)

pchic.to_csv(args.hg38_PCHiC, sep='\t')