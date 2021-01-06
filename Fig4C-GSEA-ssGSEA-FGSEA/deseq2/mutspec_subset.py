import pandas as pd

GFP_659 = pd.read_csv("deseq2_outputs/GFP_vs_659_results_deseq2.csv",index_col=0)
WT_659 = pd.read_csv("deseq2_outputs/WT_vs_659_results_deseq2.csv",index_col=0)

mutspec_genes =  pd.read_csv("combined_kallisto_abundance_genes_filtered_transformed_RNF43_G659fs_mutspec.tsv",sep="\t",index_col="#gene_id").index.tolist()

print(GFP_659.head())
print(WT_659.head())

print(len(mutspec_genes))

GFP_659[GFP_659.index.isin(mutspec_genes)].to_csv("GFP_vs_659_results_deseq2_mutspec.csv")
WT_659[WT_659.index.isin(mutspec_genes)].to_csv("WT_vs_659_results_deseq2_mutspec.csv")
