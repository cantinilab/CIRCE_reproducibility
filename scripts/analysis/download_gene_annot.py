import pybiomart as pbm
import numpy as np
import argparse 


organisms = {
    "hsapiens": snakemake.output["human_annot_file"],
    "mmusculus": snakemake.output["mouse_annot_file"]
}

for organism in organisms:
    if organism == 'hsapiens':
        dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')
    if organism == 'mmusculus':
        dataset = pbm.Dataset(name='mmusculus_gene_ensembl',  host='http://www.ensembl.org')

    annot = dataset.query(
        attributes=[
            'chromosome_name',               # Chromosome
            'transcription_start_site',     # TSS
            'start_position',               # Gene body start
            'end_position',                 # Gene body end
            'strand',                       # Strand
            'external_gene_name',           # Gene name
            'transcript_biotype'            # Biotype
        ]
    )
    annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].to_numpy(dtype = str)
    filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
    annot.columns = ['Chromosome', 'TSS', 'Start', 'End', 'Strand', 'Gene', 'Transcript_type']
    annot = annot[annot.Transcript_type == 'protein_coding']
    annot["Strand"] = annot["Strand"].replace({1: "+", -1: "-"})
    annot.Start = annot.Start.astype(np.int32)
    annot.End = annot.End.astype(np.int32)
    annot.TSS = annot.TSS.astype(np.int32)

    annot.to_csv(
        organisms[organism], sep="\t", index=False)