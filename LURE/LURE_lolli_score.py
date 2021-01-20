import os
import itertools
import seaborn as sns
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from lifelines.statistics import logrank_test
from lifelines.statistics import survival_difference_at_fixed_point_in_time_test



lure = pd.read_csv("outputs/giannakis_updated_LURE_matrix.csv",
                    index_col = 0).transpose()
lure.index = [i.replace(".","-") for i in lure.index]


RNF43_catches = lure[lure["RNF43_TRUNCATING"]=="TRUNCATING"].index.tolist()
RNF43_catches = [i.replace(".","-") for i in RNF43_catches]
print("Total number of RNF43 catches: ", len(RNF43_catches))

BRAF_catches = lure[lure["BRAF_MISSENSE"]=="MISSENSE"].index.tolist()
BRAF_catches = [i.replace(".","-") for i in BRAF_catches]

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

##########################
# LURE classifier scores
##########################

lure_score = pd.read_csv("inputs/V60_TCGA_COAD_BRAF_MISSENSE_set_cover_classifier_scores_2019_06_24_12_37_mutation_specific_cytoband_fusion_5_15_2019.csv",
                        index_col = 0)

lure_score.index = [i.replace(".","-") for i in lure_score.index]
score_mut = cbio_giannakis.join(lure_score, how='outer') # add LURE score
score_mut = score_mut.join(lure, how='outer') # add LURE matrix

# print(score_mut.sort_values(by="x").to_string())

#which samples are not catches?
print("Not catches (low LURE score):")
print(score_mut[(score_mut["x"]<0.5) & (score_mut.index.isin(RNF43_catches))])

#############################################
#  scores in RNF43 events without BRAF
############################################
print("\n")
print(" ------------- scores in RNF43 events without BRAF:")

#get samples that with RNF43 events but not BRAF events
score_mut_RNF43_noBRAF =score_mut[(score_mut["RNF43_TRUNCATING"]=="TRUNCATING") & (score_mut["BRAF_MISSENSE"]!="MISSENSE")]
# score_mut =score_mut[(score_mut["RNF43_TRUNCATING"]=="TRUNCATING") ]

print(score_mut_RNF43_noBRAF)

print("\n")
print("only G659fs")
print(score_mut_RNF43_noBRAF[score_mut_RNF43_noBRAF["RNF43"].str.join(sep=",")=="G659fs"]["x"].shape[0])
print(score_mut_RNF43_noBRAF[score_mut_RNF43_noBRAF["RNF43"].str.join(sep=",")=="G659fs"]["x"].mean())

print("\n")
print("only R117fs")
print(score_mut_RNF43_noBRAF[score_mut_RNF43_noBRAF["RNF43"].str.join(sep=",")=="R117fs"]["x"].shape[0])
print(score_mut_RNF43_noBRAF[score_mut_RNF43_noBRAF["RNF43"].str.join(sep=",")=="R117fs"]["x"].mean())

print("\n")
print("both R117fs and G659fs")
print(score_mut_RNF43_noBRAF[score_mut_RNF43_noBRAF["RNF43"].apply(str).str.contains("R117fs","G659fs")]["x"].shape[0])
print(score_mut_RNF43_noBRAF[score_mut_RNF43_noBRAF["RNF43"].apply(str).str.contains("R117fs","G659fs")]["x"].mean())

print("\n")
print("not R117fs or G659fs:")
print(score_mut_RNF43_noBRAF[~score_mut_RNF43_noBRAF["RNF43"].apply(str).str.contains("R117fs","G659fs")]["x"].shape[0])
print(score_mut_RNF43_noBRAF[~score_mut_RNF43_noBRAF["RNF43"].apply(str).str.contains("R117fs","G659fs")]["x"].mean())

#############################################
#  scores in SINGLE - RNF43 events without BRAF
############################################
print("\n")
print(" ------------- scores in SINGLE - RNF43 events without BRAF:")
score_mut_RNF43_noBRAF_unique = score_mut_RNF43_noBRAF[score_mut_RNF43_noBRAF["RNF43"].str.len()==1]
print(score_mut_RNF43_noBRAF_unique)

print("\n")
print("only G659fs")
print(score_mut_RNF43_noBRAF_unique[score_mut_RNF43_noBRAF_unique["RNF43"].str.join(sep=",")=="G659fs"]["x"].shape[0])
print(score_mut_RNF43_noBRAF_unique[score_mut_RNF43_noBRAF_unique["RNF43"].str.join(sep=",")=="G659fs"]["x"].mean())

print("\n")
print("only R117fs")
print(score_mut_RNF43_noBRAF_unique[score_mut_RNF43_noBRAF_unique["RNF43"].str.join(sep=",")=="R117fs"]["x"].shape[0])
print(score_mut_RNF43_noBRAF_unique[score_mut_RNF43_noBRAF_unique["RNF43"].str.join(sep=",")=="R117fs"]["x"].mean())

print("\n")
print("not R117fs or G659fs:")
print(score_mut_RNF43_noBRAF_unique[~score_mut_RNF43_noBRAF_unique["RNF43"].apply(str).str.contains("R117fs","G659fs")]["x"].shape[0])
print(score_mut_RNF43_noBRAF_unique[~score_mut_RNF43_noBRAF_unique["RNF43"].apply(str).str.contains("R117fs","G659fs")]["x"].mean())


################################################
# plot classification score distribution
################################################

# score_mut_RNF43_noBRAF_unique["RNF43"]= score_mut_RNF43_noBRAF_unique["RNF43"].apply(str)
score_mut_RNF43_noBRAF_unique["RNF43"]= score_mut_RNF43_noBRAF_unique["RNF43"].str.join(sep=",")

plt.figure(figsize=(10,1.5))
ax = plt.subplot(111)
ax = sns.stripplot(x="RNF43", y="x", data=score_mut_RNF43_noBRAF_unique,
                    color="black", jitter=False, s=12, alpha=.8,
                    order = [ 'L17*',  'R117fs', 'R145*', 'A629fs', 'G659fs'])

plt.axhline(y=0.5, color = "red")
# plt.xticks(rotation=90)
plt.ylim((.3,.85))

plt.xlabel('', fontsize=18)
plt.ylabel('Classifier Score', fontsize=16)

plt.savefig("outputs/stripplot_lure_score_single_RNF43_noBRAF.png",dpi=300,bbox_inches="tight")
plt.savefig("outputs/stripplot_lure_score_single_RNF43_noBRAF.svg",dpi=300,bbox_inches="tight")

plt.clf()
plt.cla()


##################################################
# lollipop for all RNF43 samples LURED (over 0.5)
##################################################
RNF43_catches = score_mut[(score_mut["x"]>0.5) & (score_mut.index.isin(RNF43_catches))]["RNF43"].tolist()
#flatten list of all RNF43 mutations in LURED samples
RNF43_catches_flat =  [item for sublist in RNF43_catches for item in sublist]

# print("Total # of RNF43 truncating variants:", len(lured_rnf43_muts_flat))
print("Total # of RNF43 truncating variants:", len(RNF43_catches_flat))
count  = Counter(RNF43_catches_flat)

for key, value in count.items():
     print(key, value)

cmd = "../Fig1A-RNF43-lollipops/lollipops-v1.5.2-mac64/lollipops  -domain-labels=off  -w=700 -o outputs/lolli_all_sig.svg RNF43 "

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

##################################################
# lollipop for all RNF43 samples LURED (all)
##################################################
RNF43_catches = score_mut[score_mut.index.isin(RNF43_catches)]["RNF43"].tolist()
#flatten list of all RNF43 mutations in LURED samples
RNF43_catches_flat =  [item for sublist in RNF43_catches for item in sublist]

# print("Total # of RNF43 truncating variants:", len(lured_rnf43_muts_flat))
print("Total # of RNF43 truncating variants:", len(RNF43_catches_flat))
count  = Counter(RNF43_catches_flat)

for key, value in count.items():
     print(key, value)

cmd = "../Fig1A-RNF43-lollipops/lollipops-v1.5.2-mac64/lollipops  -domain-labels=off  -w=700 -o outputs/lolli_all.svg RNF43 "

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

############################################################
# lollipop for RNF43 unique and no BRAF samples (over 0.5)
#############################################################
score_mut_RNF43_noBRAF_unique_catches = score_mut[(score_mut["x"]>0.5) &
                                        (score_mut["RNF43_TRUNCATING"]=="TRUNCATING") &
                                        (score_mut["BRAF_MISSENSE"]!="MISSENSE") &
                                        (score_mut["RNF43"].str.len()==1) ]["RNF43"].tolist()



RNF43_noBRAF_unique_catches_flat =  [item for sublist in score_mut_RNF43_noBRAF_unique_catches for item in sublist]

count  = Counter(RNF43_noBRAF_unique_catches_flat)

for key, value in count.items():
     print(key, value)

cmd = "../Fig1A-RNF43-lollipops/lollipops-v1.5.2-mac64/lollipops  -domain-labels=off  -w=700 -o outputs/lolli_noBRAF_unique_sig.svg RNF43 "

#plot only * truncating
for key, val in count.items():
    if val == 1:
        size=1
    else:
        size=val**2
    color = "#000000"
    cmd += str(key)+color+"@"+str(size)+" "

print(cmd)
os.system(cmd)

############################################################
# lollipop for RNF43 unique and no BRAF samples (all)
#############################################################
score_mut_RNF43_noBRAF_unique_catches = score_mut[
                                        (score_mut["RNF43_TRUNCATING"]=="TRUNCATING") &
                                        (score_mut["BRAF_MISSENSE"]!="MISSENSE") &
                                        (score_mut["RNF43"].str.len()==1) ]["RNF43"].tolist()



RNF43_noBRAF_unique_catches_flat =  [item for sublist in score_mut_RNF43_noBRAF_unique_catches for item in sublist]

count  = Counter(RNF43_noBRAF_unique_catches_flat)

for key, value in count.items():
     print(key, value)

cmd = "../Fig1A-RNF43-lollipops/lollipops-v1.5.2-mac64/lollipops  -domain-labels=off  -w=700 -o outputs/lolli_noBRAF_unique_all.svg RNF43 "

#plot only * truncating
for key, val in count.items():
    if val == 1:
        size=1
    else:
        size=val**2
    color = "#000000"
    cmd += str(key)+color+"@"+str(size)+" "

print(cmd)
os.system(cmd)
