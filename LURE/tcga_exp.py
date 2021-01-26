import os
import itertools
import seaborn as sns
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

#from LURE_lolli_
score_mut = pd.read_csv("outputs/LUREscore_mut_status.csv", index_col=0)

#Xena - gene-leveltranscription estimates, as in log2(x+1) transformed RSEM normalized count
tcga_exp = pd.read_csv("inputs/HiSeqV2", sep="\t",index_col=0)
tcga_exp.columns =  [i[:-3] for i in tcga_exp.columns]

"""
#import kegg gene set
mapk_genes = None
wnt_genes = None
cancer_genes = None
cytokine_genes = None
with open("inputs/c2.cp.kegg.v7.2.symbols.gmt", "r") as kegg:
    for line in kegg:
        if line.split()[0] =="KEGG_MAPK_SIGNALING_PATHWAY":
            mapk_genes = line.split()[2:]
        if line.split()[0] =="KEGG_WNT_SIGNALING_PATHWAY":
            wnt_genes = line.split()[2:]
        if line.split()[0] =="KEGG_PATHWAYS_IN_CANCER":
            cancer_genes = line.split()[2:]
        if line.split()[0] =="KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION":
            cytokine_genes = line.split()[2:]

genes = []
with open("inputs/h.all.v6.0.symbols.gmt", "r") as hm:
    for line in hm:
        # if line.split()[0] in ["HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_HYPOXIA","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]:

        if line.split()[0] in ["HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN"]:
            genes.extend(line.split()[2:])
"""

#https://www.nature.com/articles/s41698-018-0051-4 "MAPK Pathway Activity Score "
tengenes = "SPRY2, SPRY4, ETV4, ETV5, DUSP4, DUSP6, CCND1, EPHA2,  EPHA4".replace(',','').split()
wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3"]
# ssh= ["GLI1","GLI2"]
#EMT from https://bmccancer.biomedcentral.com/articles/10.1186/s12885-020-07615-5/figures/5
emt = "PLS3 TRAP1 IGFBP3 CLU FOSL1 FSTL3 SNAI1 CCL19 HAND1 PCOLCE2 SCG2".split()

# hypoxiia https://cancerci.biomedcentral.com/articles/10.1186/s12935-019-0964-1
hypoxia = "ORAI3 FBP1 VEGFA MVD GAS6 CYB5R3 ZBTB44 MDM2 PRELID2 CCNG1 FAM117B CASP6 TRAF3 PRP1B".split()

# tcga_exp = tcga_exp[tcga_exp.index.isin(mapk_genes+wnt_genes+cancer_genes+cytokine_genes+genes)]
# tcga_exp = tcga_exp[tcga_exp.index.isin(mapk_genes+genes+tengenes)]
# tcga_exp = tcga_exp[tcga_exp.index.isin(tengenes+wnt+ssh)]
tcga_exp = tcga_exp[tcga_exp.index.isin(tengenes+wnt)]


#remove low std genes
# tcga_exp = tcga_exp.loc[tcga_exp.std(axis=1) > 1, :]

print(tcga_exp)

# sys.exit()
# print(score_mut["RNF43_TRUNCATING"].transpose().tolist())

tcga_exp = tcga_exp.loc[:,~tcga_exp.columns.duplicated()]
print("drop dup")
print(tcga_exp)

#find tcga samples in both
samples = [i for i in tcga_exp.columns if i in score_mut.index ]
tcga_exp = tcga_exp[samples]
score_mut = score_mut[score_mut.index.isin(samples)][["RNF43_TRUNCATING","BRAF_MISSENSE"]]
tcga_exp = tcga_exp.append(score_mut[["RNF43_TRUNCATING","BRAF_MISSENSE"]].transpose())

rnf43_status = tcga_exp.loc["RNF43_TRUNCATING"].tolist()
braf_status = tcga_exp.loc["BRAF_MISSENSE"].tolist()




#convert mut status to color
rnf43_status_c = []
for i in rnf43_status:
    if i == "TRUNCATING":
        rnf43_status_c.append("red")
    else:
        rnf43_status_c.append("grey")
# print(rnf43_status_c)
# print(len(rnf43_status_c))

#add color for left right side
#add color for BRAF score


#convert mut status to
braf_status_c = []
for i in braf_status:
    if i == "MISSENSE":
        braf_status_c.append("blue")
    else:
        braf_status_c.append("grey")
# print(braf_status_c)
# print(len(braf_status_c))

tcga_exp = tcga_exp.drop("RNF43_TRUNCATING", axis=0)
tcga_exp = tcga_exp.drop("BRAF_MISSENSE", axis=0)
tcga_exp = tcga_exp.convert_objects(convert_numeric=True)

# for i in ["single","complete", "average", "weighted", "centroid","median","ward"]:

print(tcga_exp.index.tolist())

gene_colors = []
for i in tcga_exp.index.tolist():
    print(i)
    if i in tengenes:
        gene_colors.append("lightblue")
    elif i in wnt:
        gene_colors.append("pink")
    else:
        gene_colors.append("grey")




sns.clustermap(tcga_exp,z_score=0,cmap="coolwarm",xticklabels=False,
            yticklabels=True,col_colors=[rnf43_status_c,braf_status_c],
            row_colors=gene_colors, method = "ward")




handles = [Patch(facecolor=lut[name]) for name in lut]
plt.legend(handles, lut, title='Genes',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')


plt.tight_layout()
plt.savefig("outputs/heatmap.png",bbox_inches = "tight",dpi=400)
plt.clf()
plt.close()
