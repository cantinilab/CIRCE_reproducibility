import pandas as pd

keys = ["circe_singlecell", "cicero_singlecell", "cicero_pseudocell", "circe_pseudocell"]
paths = [snakemake.input[key] for key in keys]

l_df = [pd.read_csv(path, sep="\t") for path in paths]

for i in range(len(l_df)):
    l_df[i] = l_df[i][l_df[i]["coaccess"]!=0]
    l_df[i] = l_df[i].rename(columns={"coaccess": keys[i]},).set_index(["Peak1", "Peak2"])
    l_df[i] = l_df[i][keys[i]]

df_all = pd.concat(l_df, axis=1)
df_all.corr(method="pearson").to_csv(snakemake.output[0])
df_all.corr(method="spearman").to_csv(snakemake.output[1])
