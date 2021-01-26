import os
import itertools
import seaborn as sns
import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

#from LURE_lolli_
score_mut = pd.read_csv("outputs/LUREscore_mut_status.csv", index_col=0)

#Xena - gene-leveltranscription estimates, as in log2(x+1) transformed RSEM normalized count
# tcga_exp = pd.read_csv("inputs/HiSeqV2", sep="\t",index_col=0)
tcga_exp = pd.read_csv("inputs/XenaCOADREAD_HiSeqV2", sep="\t",index_col=0)
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
mapkgenes = "SPRY2, SPRY4, ETV4, ETV5, DUSP4, DUSP6, CCND1, EPHA2,  EPHA4".replace(',','').split()
wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3"]
# ssh= ["GLI1","GLI2"]
#EMT from https://bmccancer.biomedcentral.com/articles/10.1186/s12885-020-07615-5/figures/5
emt = "PLS3 TRAP1 IGFBP3 CLU FOSL1 FSTL3 SNAI1 CCL19 HAND1 PCOLCE2 SCG2".split()

# hypoxiia https://cancerci.biomedcentral.com/articles/10.1186/s12935-019-0964-1
hypoxia = "ORAI3 FBP1 VEGFA MVD GAS6 CYB5R3 ZBTB44 MDM2 PRELID2 CCNG1 FAM117B CASP6 TRAF3 PRP1B".split()

tcga_exp = tcga_exp[tcga_exp.index.isin(mapkgenes+wnt)]

#remove low std genes
# tcga_exp = tcga_exp.loc[tcga_exp.std(axis=1) > 1, :]

#drop duplicated samples
tcga_exp = tcga_exp.loc[:,~tcga_exp.columns.duplicated(keep="first")]


#find tcga samples in both
samples = [i for i in tcga_exp.columns if i in score_mut.index ]
tcga_exp = tcga_exp[samples]
score_mut = score_mut[score_mut.index.isin(samples)][["RNF43_TRUNCATING","BRAF_MISSENSE","x","RNF43"]]
tcga_exp = tcga_exp.append(score_mut[["RNF43_TRUNCATING","BRAF_MISSENSE","x","RNF43"]].transpose())

##################################
# adding color columns
##################################

rnf43_status = tcga_exp.loc["RNF43_TRUNCATING"].tolist()
braf_status = tcga_exp.loc["BRAF_MISSENSE"].tolist()
lure_score = tcga_exp.loc["x"].tolist()

#convert RNF43 mut status to color
rnf43_status_c = []
for i in rnf43_status:
    if i == "TRUNCATING":
        rnf43_status_c.append("orange")
    else:
        rnf43_status_c.append("#f0edf2")

#convert braf mut status to color
braf_status_c = []
for i in braf_status:
    if i == "MISSENSE":
        braf_status_c.append("green")
    else:
        braf_status_c.append("#f0edf2")


#adding scale for LURE classifier score
norm = matplotlib.colors.Normalize(vmin=0, vmax=1, clip=True)
mapper = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.YlGn)
tcga_exp.loc['x_color'] = tcga_exp.loc['x'].apply(lambda x: mcolors.to_hex(mapper.to_rgba(x)))
lure_scale_c = tcga_exp.loc['x_color'].tolist()


#left vs right side
clinical = pd.read_csv("inputs/coadread_tcga_clinical_data.tsv", sep="\t",
                        index_col = "Patient ID")

left = ['Splenic Flexure','Descending Colon','Sigmoid Colon','Rectosigmoid Junction','Rectum']
right = ['Cecum','Ascending Colon','Hepatic Flexure','Transverse Colon']

left_samples = clinical[clinical["Patient Primary Tumor Site"].isin(left)].index.tolist()
right_samples = clinical[clinical["Patient Primary Tumor Site"].isin(right)].index.tolist()

#convert side  status to color
side_status_c = []
for i in tcga_exp.columns:
    if i in  left_samples:
        side_status_c.append("lightblue")
    elif i in right_samples :
        side_status_c.append("#804ea3")
    else:
        side_status_c.append("#f0edf2")


"""
#clinical from Xena that has MSI status
xenaclinical = pd.read_csv("inputs/COADREAD_clinicalMatrix", sep="\t",
                        index_col = "_PATIENT")
xenaclinical = xenaclinical[xenaclinical.index.isin(tcga_exp.columns)]["MSI_updated_Oct62011"]

print(xenaclinical.value_counts())
"""

#prep df for clustermap
tcga_exp = tcga_exp.drop("RNF43_TRUNCATING", axis=0)
tcga_exp = tcga_exp.drop("BRAF_MISSENSE", axis=0)
tcga_exp = tcga_exp.drop("x", axis=0)
tcga_exp = tcga_exp.drop("x_color", axis=0)
tcga_exp = tcga_exp.drop("RNF43", axis=0)
tcga_exp = tcga_exp.convert_objects(convert_numeric=True)


#now with only genes rows left, label colors
gene_colors = []
for i in tcga_exp.index.tolist():
    if i in mapkgenes:
        gene_colors.append("lightgrey")
    elif i in wnt:
        gene_colors.append("darkgrey")
    else:
        gene_colors.append("pink")


# for i in ["single","complete", "average", "weighted", "centroid","median","ward"]:
sns.clustermap(tcga_exp,z_score=0,cmap="coolwarm",xticklabels=False,
            yticklabels=True,col_colors=[lure_scale_c,braf_status_c,rnf43_status_c],
            row_colors=gene_colors, method = "ward", colors_ratio = .03, figsize = (12,8))


# handles = [Patch(facecolor=lut[name]) for name in lut]
# plt.legend(handles, lut, title='Genes',
#            bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')


plt.tight_layout()
plt.savefig("outputs/heatmap.png",bbox_inches = "tight",dpi=400)
plt.clf()
plt.close()
