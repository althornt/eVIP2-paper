import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from collections import Counter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

import os
import itertools

lure = pd.read_csv("outputs/giannakis_updated_LURE_matrix.csv",
                    index_col = 0).transpose()

RNF43_catches = lure[lure["RNF43_TRUNCATING"]=="TRUNCATING"].index.tolist()
RNF43_catches = [i.replace(".","-") for i in RNF43_catches]
print("Total number of RNF43 catches: ", len(RNF43_catches))


# cbioportal samples RNF43 status
RNF43_muts = pd.read_csv("inputs/RNF43-mutations-cbioportal.txt", sep="\t",
                        index_col = "SAMPLE_ID")
RNF43_muts.index = [i[:-3] for i in RNF43_muts.index]

#lured samples
lured_rnf43_muts = RNF43_muts[RNF43_muts.index.isin(RNF43_catches)]["RNF43"].tolist()
cbioportal = RNF43_muts[RNF43_muts.index.isin(RNF43_catches)]["RNF43"].to_frame()

cbioportal["RNF43"] = cbioportal["RNF43"].str.split()

#remove WT at and non truncating variants from cbioportal
for index, row in cbioportal.iterrows():
    keep = []
    for i in row["RNF43"]:
        if i != "WT" and "*"  in i :
            #make mut nomenclature similar to giannakis
            if "fs*"  in i :
                keep.append(i.split("fs*")[0][:-1]+"fs")
            else:
                keep.append(i)
    row["RNF43"] = keep

# giannakis RNF43 calls
giannakis = pd.read_csv("inputs/Giannakis-TCGA-COAD-RNF43.txt", sep="\t",index_col = "Tumor_Sample_Barcode")
giannakis.index = [i[-12:] for i in giannakis.index]
# giannakis.index = [i.replace("-",".") for i in giannakis.index]
giannakis_fs = giannakis[(giannakis["Protein_Change"].str.contains("fs")) | (giannakis["Protein_Change"].str.endswith("*"))]
giannakis_fs["Protein_Change"] = [str(x).replace('p.','') for x in giannakis_fs["Protein_Change"]]

#combine calls from each sample
giannakis_fs_gp = giannakis_fs['Protein_Change'].groupby(giannakis_fs.index).apply(list).to_frame()

#combine cbioportal and giannakis df
cbio_giannakis = cbioportal.join(giannakis_fs_gp["Protein_Change"], how='outer')

#reduce to LURE caatches
cbio_giannakis = cbio_giannakis[cbio_giannakis.index.isin(RNF43_catches)]

#combining columns
for index, row in cbio_giannakis.iterrows():
    try:
        row['RNF43'].extend(row['Protein_Change'])
        row['RNF43'] = list(set(row['RNF43']))

    except: #only when trying to add nan to list
        pass


cbio_giannakis["RNF43"].to_csv("outputs/giannakais_cbioportal_RNF43_truncating_variants.csv")

print(cbio_giannakis["RNF43"])


lured_RNF43_all_counts = cbio_giannakis["RNF43"].tolist()


#######################
# lollipop
#######################

#flatten list of all RNF43 mutations in LURED samples
lured_rnf43_muts_flat =  [item for sublist in lured_RNF43_all_counts for item in sublist]

print("Total # of RNF43 truncating variants:", len(lured_rnf43_muts_flat))

count  = Counter(lured_rnf43_muts_flat)

for key, value in count.items():
     print(key, value)



cmd = "../Fig1A-RNF43-lollipops/lollipops-v1.5.2-mac64/lollipops  -domain-labels=off  -w=700 -o outputs/lolli_out_update.svg RNF43 "

#plot only * truncating
for key, val in count.items():
    if val == 1:
        size=1
    else:
        size=val**4


    color = "#000000"

    cmd += str(key)+color+"@"+str(size)+" "

print(cmd)
os.system(cmd)
