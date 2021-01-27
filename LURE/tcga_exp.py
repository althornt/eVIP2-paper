import os
import itertools
import seaborn as sns
import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


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
mapkgenes = "PHLDA1, SPRY2, SPRY4, ETV4, ETV5, DUSP4, DUSP6, CCND1, EPHA2,  EPHA4".replace(',','').split()
# wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3"]
# wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3","DKK1","WIF1","LRP6","CDH1","CTNNB1","GSK3B"]
wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3","WIF1"]

# wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3","DKK1","WIF1"]

mmr = ["MLH1"]

"""
# wnt2 = ["DKK1","APCDD1","WIF1","ZNRF3"]
# wnt2 = ["DKK1","WIF1","LRP6"]

# ssh= ["GLI1","GLI2"]
#EMT from https://bmccancer.biomedcentral.com/articles/10.1186/s12885-020-07615-5/figures/5
emt = "PLS3 TRAP1 IGFBP3 CLU FOSL1 FSTL3 SNAI1 CCL19 HAND1 PCOLCE2 SCG2".split()

# hypoxiia https://cancerci.biomedcentral.com/articles/10.1186/s12935-019-0964-1
hypoxia = "ORAI3 FBP1 VEGFA MVD GAS6 CYB5R3 ZBTB44 MDM2 PRELID2 CCNG1 FAM117B CASP6 TRAF3 PRP1B".split()
"""
tcga_exp = tcga_exp[tcga_exp.index.isin(mapkgenes+wnt+mmr)]



#remove low std genes
# tcga_exp = tcga_exp.loc[tcga_exp.std(axis=1) > 1, :]

#drop duplicated samples
tcga_exp = tcga_exp.loc[:,~tcga_exp.columns.duplicated(keep="first")]


#find tcga samples in both
samples = [i for i in tcga_exp.columns if i in score_mut.index ]
tcga_exp = tcga_exp[samples]
score_mut = score_mut[score_mut.index.isin(samples)][["RNF43_TRUNCATING","BRAF_MISSENSE","x","RNF43"]]
#rename RNF43 so its not confused with gene expression row
score_mut = score_mut.rename(columns={"RNF43": "RNF43_muts"})
tcga_exp = tcga_exp.append(score_mut[["RNF43_TRUNCATING","BRAF_MISSENSE","x","RNF43_muts"]].transpose())

##################################
# adding color columns
##################################

rnf43_status = tcga_exp.loc["RNF43_TRUNCATING"].tolist()
braf_status = tcga_exp.loc["BRAF_MISSENSE"].tolist()
lure_score = tcga_exp.loc["x"].tolist()
RNF43_muts  = tcga_exp.loc["RNF43_muts"].tolist()

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
        braf_status_c.append("darkgreen")
    else:
        braf_status_c.append("#f0edf2")

RNF43_mut_type_c = []
for i in RNF43_muts:
    if "G659fs" in str(i):
        RNF43_mut_type_c.append("gold")
    else:
        RNF43_mut_type_c.append("#f0edf2")

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
tcga_exp = tcga_exp.drop("RNF43_muts", axis=0)
tcga_exp = tcga_exp.convert_objects(convert_numeric=True)


#now with only genes rows left, label colors
gene_colors = []
for i in tcga_exp.index.tolist():
    if i in mapkgenes:
        gene_colors.append("slateblue")
    elif i in wnt:
        gene_colors.append("slategrey")
    elif i in mmr:
        gene_colors.append("plum")

    else:
        gene_colors.append("pink")



# for i in ["single","complete", "average", "weighted", "centroid","median","ward"]:
# sns.clustermap(tcga_exp,z_score=0,center=0,cmap="coolwarm",xticklabels=False,
#             yticklabels=True,row_colors=gene_colors, method = "ward",
#             colors_ratio = .02, figsize = (12,8),
#             col_colors=[lure_scale_c,braf_status_c,rnf43_status_c,RNF43_mut_type_c],
#             row_cluster=True, cbar_pos=None,vmin=-5, vmax=4)

sns.clustermap(tcga_exp,z_score=0,center=0,cmap="bwr",xticklabels=False,
            yticklabels=True,row_colors=gene_colors, method = "ward",
            colors_ratio = .02, figsize = (12,8),
            col_colors=[lure_scale_c,braf_status_c,rnf43_status_c,RNF43_mut_type_c],
            row_cluster=True)

plt.tight_layout()
plt.savefig("outputs/tcga_heatmap.png",bbox_inches = "tight",dpi=400)
plt.savefig("outputs/tcga_heatmap.svg",bbox_inches = "tight",dpi=400)

plt.clf()
plt.close()

##################################################################
# Putting legends and colorbars in another file
##################################################################

mapk_patch = mpatches.Patch(color='slateblue', label='MAPK genes')
wnt_patch = mpatches.Patch(color='slategrey', label='WNT genes')
mmr_patch = mpatches.Patch(color='plum', label='MMR gene')
plt.legend(handles=[wnt_patch, mapk_patch, mmr_patch],loc='upper center')

#color baar for LURE classifier score
sm = plt.cm.ScalarMappable(cmap=plt.cm.YlGn, norm=plt.Normalize(vmin=0, vmax=1))
sm._A = []
plt.colorbar(sm)

#color baar for heatmap zscore vaalue
sm = plt.cm.ScalarMappable(cmap=plt.cm.bwr, norm=plt.Normalize(vmin=-5, vmax=5))
sm._A = []
plt.colorbar(sm)

plt.savefig("outputs/tcga_heatmap_legend.svg",bbox_inches = "tight",dpi=400)
plt.clf()
plt.close()
