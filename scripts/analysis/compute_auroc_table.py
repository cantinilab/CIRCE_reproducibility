from sklearn.metrics import roc_auc_score
import pandas as pd

# Compute auroc score per output and save it as a dictionary
# open dataframe
auroc_files = snakemake.input['auroc_files']
auroc_dict = {}
for auroc_file in auroc_files:
    auroc_df = pd.read_csv(auroc_file, sep='\t')
    print(auroc_file)
    auroc = roc_auc_score(auroc_df['score'], auroc_df['name'])
    auroc_dict[auroc_file] = auroc

# Save dictionary as a table
auroc_table = pd.DataFrame(auroc_dict.items(), columns=['file', 'auroc'])
auroc_table.to_csv(snakemake.output['auroc_table'], sep='\t', index=False)