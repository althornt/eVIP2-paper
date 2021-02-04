import pandas as pd
import matplotlib
import seaborn as sns
import sys
import matplotlib.pyplot as plt
import numpy as  np
from scipy import stats



#cbpioportal RSEM
tcga_exp = pd.read_csv("inputs/cbioportal_data_RNA_Seq_v2_expression_median.txt", sep="\t",index_col=0)
tcga_exp =  tcga_exp.drop(columns="Entrez_Gene_Id")
tcga_exp.columns =  [i[:-3] for i in tcga_exp.columns]
tcga_exp = np.log2(tcga_exp+1) #log2transform
#z-score expression data
tcga_exp_ = stats.zscore(tcga_exp, axis=1)
tcga_exp = pd.DataFrame(tcga_exp_,
            columns = tcga_exp.columns,
            index = tcga_exp.index)


#BRAF-scores and RNF43 mutants generated from LURE lolli script
score_mut = pd.read_csv("outputs/LUREscore_mut_status.csv", index_col=0)

#adding driver status of all mapk/rtk genes
mapk_driver_status = pd.read_csv("inputs/cbioportal-MAPK-RTK-genes-driver-sample-matrix.txt",
                sep="\t", index_col = 0)
mapk_driver_status.index =  [i.split(":")[1][:-3] for i in mapk_driver_status.index]

#samples in all three df
samples =  set(score_mut.index ).intersection(tcga_exp.columns, mapk_driver_status.index)
score_mut = score_mut[score_mut.index.isin(samples)][["RNF43_TRUNCATING","BRAF_MISSENSE","x"]]
mapk_driver_status = mapk_driver_status[mapk_driver_status.index.isin(score_mut.index)]
score_mut = score_mut.merge(mapk_driver_status, how='outer', left_index=True, right_index=True)
score_mut = score_mut.drop(columns=["Altered","x"])


RNF43_WT_BRAF_WT = score_mut[ (score_mut["RNF43_TRUNCATING"] != "TRUNCATING")
            & (score_mut["BRAF_MISSENSE"] != "MISSENSE") ].drop(columns=["RNF43_TRUNCATING","BRAF_MISSENSE"])

print("RNF43_WT_BRAF_WT")
# print(RNF43_WT_BRAF_WT.loc[:,RNF43_WT_BRAF_WT.sum() !=0].sum().sort_values(ascending=False))
print("total samples:",RNF43_WT_BRAF_WT.shape[0])
print(RNF43_WT_BRAF_WT.loc[:,RNF43_WT_BRAF_WT.sum() !=0].sum().sort_values(ascending=False)/RNF43_WT_BRAF_WT.shape[0])
print(RNF43_WT_BRAF_WT.loc[RNF43_WT_BRAF_WT.sum(axis=1)<1])

#RNF43 WT & BRAF Missense
RNF43_WT_BRAF_MISS = score_mut[ (score_mut["RNF43_TRUNCATING"] != "TRUNCATING")
            & (score_mut["BRAF_MISSENSE"] == "MISSENSE")].drop(columns=["RNF43_TRUNCATING","BRAF_MISSENSE"])

print("RNF43_WT_BRAF_MISS")
# print(RNF43_WT_BRAF_MISS.loc[:,RNF43_WT_BRAF_MISS.sum() !=0].sum().sort_values(ascending=False))
print("total samples:",RNF43_WT_BRAF_MISS.shape[0])
print(RNF43_WT_BRAF_MISS.loc[:,RNF43_WT_BRAF_MISS.sum() !=0].sum().sort_values(ascending=False)/RNF43_WT_BRAF_MISS.shape[0])
print(RNF43_WT_BRAF_MISS.loc[RNF43_WT_BRAF_MISS.sum(axis=1)<1])


# RNF43 TRUNC & BRAF WT
RNF43_TRUNC_BRAF_WT = score_mut[ (score_mut["RNF43_TRUNCATING"]== "TRUNCATING")
            & (score_mut["BRAF_MISSENSE"] != "MISSENSE")].drop(columns=["RNF43_TRUNCATING","BRAF_MISSENSE"])

print("RNF43_TRUNC_BRAF_WT")
# print(RNF43_TRUNC_BRAF_WT.loc[:,RNF43_TRUNC_BRAF_WT.sum() !=0].sum().sort_values(ascending=False))
print("total samples:",RNF43_TRUNC_BRAF_WT.shape[0])
print(RNF43_TRUNC_BRAF_WT.loc[:,RNF43_TRUNC_BRAF_WT.sum() !=0].sum().sort_values(ascending=False)/RNF43_TRUNC_BRAF_WT.shape[0])
print(RNF43_TRUNC_BRAF_WT.loc[RNF43_TRUNC_BRAF_WT.sum(axis=1)<1])

# RNF43 TRUNC & BRAF MISSENSE
RNF43_TRUNC_BRAF_MISS = score_mut[(score_mut["RNF43_TRUNCATING"]== "TRUNCATING")
            & (score_mut["BRAF_MISSENSE"] == "MISSENSE")].drop(columns=["RNF43_TRUNCATING","BRAF_MISSENSE"])

print("RNF43_TRUNC_BRAF_MISS")
# print(RNF43_TRUNC_BRAF_MISS.loc[:,RNF43_TRUNC_BRAF_MISS.sum() !=0].sum().sort_values(ascending=False))
print("total samples:",RNF43_TRUNC_BRAF_MISS.shape[0])
print(RNF43_TRUNC_BRAF_MISS.loc[:,RNF43_TRUNC_BRAF_MISS.sum() !=0].sum().sort_values(ascending=False)/RNF43_TRUNC_BRAF_MISS.shape[0])
print(RNF43_TRUNC_BRAF_MISS.loc[RNF43_TRUNC_BRAF_MISS.sum(axis=1)<1])



"""
for i in [(RNF43_WT_BRAF_WT,"RNF43_WT_BRAF_WT"),(RNF43_WT_BRAF_MISS,"RNF43_WT_BRAF_MISS"),(RNF43_TRUNC_BRAF_WT,"RNF43_TRUNC_BRAF_WT"),(RNF43_TRUNC_BRAF_MISS,"RNF43_TRUNC_BRAF_MISS")]:
    df =  i[0].loc[(i[0].sum(axis=1) != 0),
                            (i[0].sum(axis=0) != 0)]

    genes_1.append(df.columns.tolist())
    sns.clustermap(df,z_score=None)
    plt.savefig("outputs/"+str(i[1])+".png",bbox_inches = "tight",dpi=400)
    plt.clf()
    plt.close()

genes_1 = [item for sublist in genes_1 for item in sublist]
print(len(set(genes_1)))
"""

#get order of mut frequency and remove genes with no driver events
score_mut_  = score_mut.drop(columns=["RNF43_TRUNCATING","BRAF_MISSENSE"])
mut_genes_by_freq = score_mut_.loc[(score_mut_.sum(axis=1) != 0),
                        (score_mut_.sum(axis=0) != 0)].sum().sort_values(ascending=False).index.tolist()

score_mut = score_mut.loc[:,score_mut.columns.isin(mut_genes_by_freq+["RNF43_TRUNCATING","BRAF_MISSENSE"])]
tcga_exp.append(score_mut.transpose())
#make order match
tcga_exp = tcga_exp[score_mut.index]


RNF43_c =  []
for i in score_mut["RNF43_TRUNCATING"].tolist():
    if i  == "TRUNCATING":
        RNF43_c.append("orange")
    else:
        RNF43_c.append("white")


#make color lists for each
l = []

l.append(RNF43_c)

for i in  mut_genes_by_freq:
    this = []
    for j in score_mut[i].tolist():
        if j == 0:
            this.append("white")
        elif j == 1:
            this.append("black")
        else:
            this.append("red")

    l.append(this)


hypoxia_genes = ["HIF1A"]
emt_genes = ["TWIST1","SNAI1"]
nfkb_genes = ["NFKB1","NFKB2","RELB"]
#https://www.nature.com/articles/s41698-018-0051-4 "MAPK Pathway Activity Score
mapkgenes = "PHLDA1, SPRY2, SPRY4, ETV4, ETV5, DUSP4, DUSP6, CCND1, EPHA2,  EPHA4".replace(',','').split()
wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3","WIF1"]
mmr = ["MLH1"]
tcga_exp = tcga_exp[tcga_exp.index.isin(mapkgenes+wnt+mmr+hypoxia_genes+nfkb_genes+emt_genes)]

"""
#full map
sns.clustermap(tcga_exp,col_colors=l,colors_ratio = .02,figsize = (12,8),xticklabels=False,cmap="bwr",
    method = "ward")
plt.savefig("outputs/testt.png",bbox_inches = "tight",dpi=400)
plt.clf()
plt.close()
"""

for i in [(RNF43_WT_BRAF_WT.index.tolist(),"RNF43_WT_BRAF_WT"),
    (RNF43_WT_BRAF_MISS.index.tolist(),"RNF43_WT_BRAF_MISS"),
    (RNF43_TRUNC_BRAF_WT.index.tolist(),"RNF43_TRUNC_BRAF_WT"),
    (RNF43_TRUNC_BRAF_MISS.index.tolist(),"RNF43_TRUNC_BRAF_MISS")]:

    #subset expression data
    df = tcga_exp[i[0]]
    #subset score mut
    mut_df = score_mut[score_mut.index.isin(i[0])]

    RNF43_c =  []
    for k in mut_df["RNF43_TRUNCATING"].tolist():
        if k  == "TRUNCATING":
            RNF43_c.append("orange")
        else:
            RNF43_c.append("white")

    BRAF_c =  []
    for p in mut_df["BRAF_MISSENSE"].tolist():
        if p == "MISSENSE":
            BRAF_c.append("green")
        else:
            BRAF_c.append("white")

    l = []
    l.append(RNF43_c)
    l.append(BRAF_c)


    mut_df = mut_df.drop(columns = ["RNF43_TRUNCATING","BRAF_MISSENSE"])
    #for each mut gene status
    for j in mut_df.columns.tolist():
        this = []
        #for each muts status in each sample
        for x in mut_df[j].tolist():

            if x == 0:
                this.append("white")
            elif x == 1:
                this.append("black")
            else:
                this.append("red")
                print(x)
        l.append(this)


    sns.clustermap(df,col_colors = l, colors_ratio = .02,figsize = (12,8),xticklabels=False,cmap="bwr",
        method = "ward")
    plt.savefig("outputs/"+str(i[1])+".png",bbox_inches = "tight",dpi=400)
    plt.clf()
    plt.close()
