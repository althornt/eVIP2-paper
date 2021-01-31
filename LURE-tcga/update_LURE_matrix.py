import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import csv



#import LURE matrix
matrix = pd.read_csv("inputs/V60_TCGA_COAD_BRAF_MISSENSE_Set_Cover_mutation_matrix_2019_06_24_12_37_mutation_specific_cytoband_fusion_5_15_2019.csv", index_col = 0).transpose()
LURE_rnf43_trunc = matrix[matrix["RNF43_TRUNCATING"]=="TRUNCATING"].index.tolist()

#getting giannakis truncating variants
giannakis = pd.read_csv("inputs/Giannakis-TCGA-COAD-RNF43.txt", sep="\t",index_col = "Tumor_Sample_Barcode")
giannakis.index = [i[-12:] for i in giannakis.index]
giannakis.index = [i.replace("-",".") for i in giannakis.index]
giannakis_fs = giannakis[(giannakis["Protein_Change"].str.contains("fs")) | (giannakis["Protein_Change"].str.endswith("*"))]
giannakis_fs = giannakis_fs["Protein_Change"].groupby(giannakis_fs.index).apply(list).to_frame()
giannakis_fs_samples = set(giannakis_fs.index.tolist())


#which giannakis calls were not used originally with LURE matrix?
#make new matrix that adds in giannakis truncating calls
new_matrix = matrix.copy()

for i in giannakis_fs_samples:
    if i not in LURE_rnf43_trunc:
        if i in matrix.index.tolist():
            new_matrix["RNF43_TRUNCATING"][i]="TRUNCATING"
        else:
            print("This TCGA sample is not in lure matrix", i)

print("RNF43 truncating calls in original matrix:", len(LURE_rnf43_trunc))
print("RNF43 truncating calls in updated matrix:",
        len(new_matrix[new_matrix["RNF43_TRUNCATING"]=="TRUNCATING"].index.tolist()))

print(matrix)
print(new_matrix)


new_matrix.transpose().to_csv("outputs/giannakis_updated_LURE_matrix.csv", quoting=csv.QUOTE_ALL)

#then run LURES oncoprint.py with the updated matrix
