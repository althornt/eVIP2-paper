import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from collections import Counter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

import os
import itertools
import seaborn as sns
import sys

lure = pd.read_csv("outputs/giannakis_updated_LURE_matrix.csv",
                    index_col = 0).transpose()

RNF43_catches = lure[lure["RNF43_TRUNCATING"]=="TRUNCATING"].index.tolist()
RNF43_catches = [i.replace(".","-") for i in RNF43_catches]
print("Total number of RNF43 catches: ", len(RNF43_catches))


# cbioportal samples RNF43 status
RNF43_muts = pd.read_csv("inputs/RNF43-mutations-cbioportal.txt", sep="\t",
                        index_col = "SAMPLE_ID")
RNF43_muts.index = [i[:-3] for i in RNF43_muts.index]

#cbioportal samples that were in lure matrix
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
# BRAF co-occurance
#######################

# cbioportal samples BRAF status
BRAF_muts = pd.read_csv("inputs/BRAF-mutations-cbioportal.txt", sep="\t",
                        index_col = "SAMPLE_ID")
BRAF_muts.index = [i[:-3] for i in BRAF_muts.index]

#LURE matrix
RNF43_catches_noBRAF = lure[(lure["RNF43_TRUNCATING"]=="TRUNCATING")
                        & (lure["BRAF_MISSENSE"].isnull())].index.tolist()
# RNF43_catches_noBRAF = [i.replace(".","-") for i in RNF43_catches_noBRAF]

print("Number of RNF43 catches without BRAF mutations:", len(RNF43_catches_noBRAF))


print(lure[lure["RNF43_TRUNCATING"]=="TRUNCATING"].shape)

lure = pd.read_csv("inputs/V60_TCGA_COAD_BRAF_MISSENSE_Set_Cover_mutation_matrix_2019_06_24_12_37_mutation_specific_cytoband_fusion_5_15_2019.csv",
                    index_col = 0).transpose()


print(lure[lure["RNF43_TRUNCATING"]=="TRUNCATING"].shape)



##########################
# LURE classifier scores
##########################

lure_score = pd.read_csv("inputs/V60_TCGA_COAD_BRAF_MISSENSE_set_cover_classifier_scores_2019_06_24_12_37_mutation_specific_cytoband_fusion_5_15_2019.csv",
                        index_col = 0)

lure_score.index = [i.replace(".","-") for i in lure_score.index]
score_mut = cbio_giannakis.join(lure_score, how='outer')


print(score_mut[(score_mut["x"]<0.5) & (score_mut.index.isin(RNF43_catches))])

score_mut["RNF43"]= score_mut["RNF43"].apply(str)


plt.figure(figsize=(4,3))
ax = plt.subplot(111)
ax = sns.boxplot(x="RNF43", y="x", data=score_mut)
plt.xticks(rotation=90)
plt.savefig("outputs/boxplot_lure_score.png",dpi=300,bbox_inches="tight")
plt.clf()
plt.cla()





def KM(KM_df, list1,label1, list2, label2):

    list1_df = KM_df[KM_df.index.isin(list1)]
    list1_df["label"] = label1

    list2_df = KM_df[KM_df.index.isin(list2)]
    list2_df["label"] = label2

    df = list1_df.append(list2_df)
    df = df.dropna()

    groups = df['label']
    ix = (groups == label1)
    T = df["Overall Survival (Months)"]
    E = df["Overall Survival Status"]
    results = logrank_test(T[~ix], T[ix], event_observed_A=E[~ix], event_observed_B=E[ix])

    p = results.p_value
    print(p)

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(4,3))
    ax = plt.subplot(111)

    for name, grouped_df in df.groupby('label'):
        print(name)
        kmf.fit(grouped_df["Overall Survival (Months)"], grouped_df["Overall Survival Status"], label=name)
        kmf.plot_survival_function(ax=ax,ci_show=False)

    ax.set_ylabel('Percent Survival')
    ax.set_xlabel('Months')


    # plt.title(tcgatype+" "+sig)
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")

    plt.text(2.,.5,"p={}".format(p), bbox={'facecolor':'w','pad':5},
         ha="right", va="top", transform=plt.gca().transAxes )
    plt.text(2.,.3,"n={}".format(len(list1)), bbox={'facecolor':'w','pad':5},
         ha="right", va="top", transform=plt.gca().transAxes )
    plt.text(2.,.1,"n={}".format(len(list2)), bbox={'facecolor':'w','pad':5},
         ha="right", va="top", transform=plt.gca().transAxes )


    plt.savefig("KM/"+label1+"_"+label2+"_KM.svg",dpi=300,bbox_inches="tight")
    plt.clf()
    plt.cla()


clinical = pd.read_csv("inputs/coadread_tcga_clinical_data.tsv", sep="\t",
                        index_col = "Patient ID")
# clinical2 = pd.read_csv("inputs/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv",
#                         sep="\t", index_col = "Patient ID")

#set up df for KM
KM_df = clinical[["Overall Survival Status","Overall Survival (Months)"]]
KM_df = KM_df.replace('0:LIVING', 0)
KM_df = KM_df.replace('1:DECEASED', 1)


#LURED RNF43 samples vs WT RNF43 samples
lured_rnf43_samples = RNF43_muts[RNF43_muts.index.isin(RNF43_catches)].index.tolist()
rnf43_wt = RNF43_muts[RNF43_muts["RNF43"]=="WT"].index.tolist()
KM(KM_df,rnf43_wt,"RNF43 WT", lured_rnf43_samples, "LURE RNF43-trunc")


#riight vs left
right = ['Cecum','Ascending Colon','Hepatic Flexure','Transverse Colon']
left = ['Splenic Flexure','Descending Colon','Sigmoid Colon','Rectosigmoid Junction','Rectum']

left_samples = clinical[clinical["Patient Primary Tumor Site"].isin(left)].index.tolist()
right_samples = clinical[clinical["Patient Primary Tumor Site"].isin(right)].index.tolist()
KM(KM_df,left_samples, "left", right_samples, "right")


# ##################################################
# # how many RNF43 samples are left vs right side
# ###############################################
#
# lured_rnf43_samples_left = [i for i in left_samples if i in lured_rnf43_samples]
# print(len(lured_rnf43_samples_left))
# lured_rnf43_samples_right = [i for i in right_samples if i in lured_rnf43_samples]
# print(len(lured_rnf43_samples_right))


braf_wt = BRAF_muts[BRAF_muts["BRAF"]=="WT"].index.tolist()
braf_v600e = BRAF_muts[BRAF_muts["BRAF"].str.contains("V600E")].index.tolist()


####################################
# specific KM comparisons
####################################

#list of samples conditions
conditions = []
for i in ["rnf43_wt","lured_rnf43_samples"]:
    for j in ["braf_wt","braf_v600e"]:
        for k in ["left_side","right_side"]:
            conditions.append([i,j,k])


dict = {"lured_rnf43_samples":lured_rnf43_samples,
        "rnf43_wt":rnf43_wt, "braf_wt":braf_wt,
        "braf_v600e":braf_v600e, "left_side":left_samples,
         "right_side":right_samples}




#combine into pairs to compare
all_km_comparisons = list(itertools.combinations(conditions,2))
for i in all_km_comparisons:
    label1= '_'.join(map(str, i[0]))
    label2= '_'.join(map(str, i[1]))

    list1 = set(dict[i[0][0]])&set(dict[i[0][1]])&set(dict[i[0][2]])
    list2 = set(dict[i[1][0]])&set(dict[i[1][1]])&set(dict[i[1][2]])


    if len(list1) > 2 and len(list2) > 2 :
        print("\n")
        print(label1, label2)
        print(len(list1),len(list2))
        KM(KM_df,list1,label1, list2, label2)


"""
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
"""
