import pybedtools
import pandas as pd
import tqdm
import anndata as ad
import circe as ci
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import os

# # Define the subfolder for temporary files
# current_dir = os.getcwd()
# temp_dir = os.path.join(current_dir, "data")
# os.environ['TMPDIR'] = temp_dir

# Open ATAC data
atac = ad.read_h5ad(snakemake.input["atac"])
atac = ci.add_region_infos(atac)
atac.var = atac.var.rename(columns = {"chromosome": "chrom"})
atac_bed = pybedtools.BedTool.from_dataframe(atac.var)

# Open PCHiC data
PCHiC = pd.read_csv(snakemake.input["PCHiC"], sep="\t")
# Rename chromosomes
PCHiC["baitChr"] = PCHiC["baitChr"] #.astype(str).apply(lambda x: f"chr{x}")
PCHiC["oeChr"] = PCHiC["oeChr"] #.astype(str).apply(lambda x: f"chr{x}")
# Ensure start and end are integers
# Keep only rows where bait and oe start are finites
# (NAs due to genomre conversion)
PCHiC = PCHiC[PCHiC["baitStart"].notna()]
PCHiC = PCHiC[PCHiC["baitEnd"].notna()]
PCHiC = PCHiC[PCHiC["oeStart"].notna()]
PCHiC = PCHiC[PCHiC["oeEnd"].notna()]
PCHiC["baitStart"] = PCHiC["baitStart"].astype(int)
PCHiC["baitEnd"] = PCHiC["baitEnd"].astype(int)
PCHiC["oeStart"] = PCHiC["oeStart"].astype(int)
PCHiC["oeEnd"] = PCHiC["oeEnd"].astype(int)

baitIDs = PCHiC["baitID"].unique()
# Add regions full names to PCHiC
PCHiC["bait_region"] = PCHiC["baitChr"] + "_" + PCHiC["baitStart"].astype(str) + "_" + PCHiC["baitEnd"].astype(str)
PCHiC["oe_region"] = PCHiC["oeChr"] + "_" + PCHiC["oeStart"].astype(str) + "_" + PCHiC["oeEnd"].astype(str)


bait_bed = pybedtools.BedTool.from_dataframe(
    PCHiC.loc[:,
        ["baitChr", "baitStart", "baitEnd"]].rename(
            columns={
                "baitChr":"chrom",
                "baitStart":"start",
                "baitEnd":"end"}
                ).drop_duplicates()
    )

# Open Circe network
circe_network = pd.read_csv(snakemake.input["predictions"], sep="\t")
circe_network["Peak1"] = circe_network["Peak1"].str.replace("-", "_")
circe_network["Peak2"] = circe_network["Peak2"].str.replace("-", "_")
circe_network = circe_network[circe_network["coaccess"] != 0]
print(circe_network.head())
circe_network_reversed = circe_network[["Peak2", "Peak1", "coaccess"]].rename(columns={"Peak2": "Peak1", "Peak1": "Peak2"})
circe_network = pd.concat([circe_network, circe_network_reversed]).reset_index(drop=True)
circe_network.columns = ["promoter", "enhancer", "coaccess"]
#circe_network = circe_network.iloc[:100000,:]
del circe_network_reversed

# Bed from all network regions
circe_regions = circe_network.loc[:, "promoter"].str.split("_", expand=True).rename(
        columns={0: "chrom", 1: "start", 2: "end"}).drop_duplicates()
circe_regions_bed = pybedtools.BedTool.from_dataframe(circe_regions)

print(circe_regions)
print(PCHiC[["baitChr", "baitStart", "baitEnd"]].drop_duplicates())

# intersect promoters with bait regions
hit_promoters = circe_regions_bed.intersect(bait_bed, wa=True, wb=True).to_dataframe().astype(str)
hit_promoters.columns = ["circe_chrom", "circe_start", "circe_end", "pchic_chrom", "pchic_start", "pchic_end"]
hit_promoters["circe_name"] = hit_promoters["circe_chrom"] + '_' + hit_promoters["circe_start"] + '_' + hit_promoters["circe_end"]
hit_promoters["pchic_name"] = hit_promoters["pchic_chrom"] + '_' + hit_promoters["pchic_start"] + '_' + hit_promoters["pchic_end"]

# filter circe_network and PCHiC to only keep links whose promoter hits:
circe_network = circe_network[circe_network["promoter"].isin(hit_promoters["circe_name"])]
circe_promoters = circe_network.loc[:, "promoter"].str.split("_", expand=True).rename(
        columns={0: "chrom", 1: "start", 2: "end"}).drop_duplicates()
PCHiC = PCHiC[PCHiC["bait_region"].isin(hit_promoters["pchic_name"])]

# Get negative links
negative_links = []
for circe_promoter in tqdm.tqdm(circe_promoters.index):
    promoter_bed = pybedtools.BedTool.from_dataframe(pd.DataFrame({
        "chrom":[circe_promoters.loc[circe_promoter, "chrom"]],
        "start":[max(0, int(circe_promoters.loc[circe_promoter, "start"])-500_000)],
        "end":[int(circe_promoters.loc[circe_promoter, "end"])+500_000]
    }))
    neg_links = atac_bed.intersect(promoter_bed).to_dataframe()
    neg_links['enhancer'] = neg_links['chrom'] + '_' + neg_links['start'].astype(str) + '_' + neg_links['end'].astype(str)
    neg_links['promoter'] = circe_promoters.loc[circe_promoter, 'chrom'] + '_' \
                            + str(circe_promoters.loc[circe_promoter, 'start']) + '_' \
                            + str(circe_promoters.loc[circe_promoter, 'end'])
    negative_links.append(neg_links[["promoter", "enhancer"]])
negative_links = pd.concat(negative_links)
circe_network = pd.concat([circe_network, negative_links]).groupby(['promoter', 'enhancer']).sum().reset_index()

# Check for each pormoter if enhancers predicted by Circe
#  is in PCHiC regions interacting with the promoter hitting region
results = []
for i in tqdm.tqdm(range(len(hit_promoters))):
    # Get Circe subnet to the promoter
    circe_subnet = circe_network[circe_network["promoter"] == hit_promoters.loc[i, "circe_name"]]
    circe_subnet_enhancers = circe_subnet["enhancer"].str.split("_", expand=True).rename(
        columns={0: "chrom", 1: "start", 2: "end"})
    circe_subnet_enhancers["coaccess"] = circe_subnet["coaccess"]
    circe_subnet_enhancers_bed = pybedtools.BedTool.from_dataframe(circe_subnet_enhancers)
    
    # Get PCHiC subnet to the promoter
    bait_subnet = PCHiC[PCHiC["bait_region"] == hit_promoters.loc[i, "pchic_name"]].loc[
        :,
        ["oeChr", "oeStart", "oeEnd"]].rename(
            columns={
                "oeChr":"chrom",
                "oeStart":"start",
                "oeEnd":"end"}
                )
    bait_subnet_bed = pybedtools.BedTool.from_dataframe(bait_subnet)
    
    # Return name and score, respectively coacces score and if the enhancer is in the PCHiC subnet
    results.append(circe_subnet_enhancers_bed.intersect(bait_subnet_bed, c=True).to_dataframe())
    results[-1]["score"] = results[-1]["score"]>0  # Binarize score

results = pd.concat(results)
results.to_csv(snakemake.output["results"], sep="\t", index=False)
