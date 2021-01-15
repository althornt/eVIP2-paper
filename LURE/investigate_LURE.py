import pandas as pd
from collections import Counter
import os

clinical = pd.read_csv("inputs/coadread_tcga_clinical_data.tsv", sep="\t",
                        index_col = "Patient ID")

clinical2 = pd.read_csv("inputs/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv",
                        sep="\t", index_col = "Patient ID")

#print cbioportal samples RNF43 status
RNF43_muts = pd.read_csv("inputs/RNF43-mutations-cbioportal.txt", sep="\t",
                        index_col = "SAMPLE_ID")
RNF43_muts.index = [i[:-3] for i in RNF43_muts.index]

#############################################
# getting RNF43 mutated samples from LURE output
#############################################

lure = pd.read_csv("inputs/V60_TCGA_COAD_BRAF_MISSENSE_Set_Cover_mutation_matrix_2019_06_24_12_37_mutation_specific_cytoband_fusion_5_15_2019.csv",
                    index_col = 0).transpose()
RNF43_catches = lure[lure["RNF43_TRUNCATING"]=="TRUNCATING"].index.tolist()
RNF43_catches = [i.replace(".","-") for i in RNF43_catches]
print("Total number of RNF43 catches: ", len(RNF43_catches))


RNF43_catches_noBRAF = lure[(lure["RNF43_TRUNCATING"]=="TRUNCATING")
                        & (lure["BRAF_MISSENSE"].isnull())].index.tolist()
RNF43_catches_noBRAF = [i.replace(".","-") for i in RNF43_catches_noBRAF]
print("Number of RNF43 catches without BRAF mutations:", len(RNF43_catches_noBRAF))
# print(RNF43_catches_noBRAF)

##############################################################
# Which RNF43 truncating mutations were caught with LURE
###########################################################

lured_rnf43_muts = RNF43_muts[RNF43_muts.index.isin(RNF43_catches)]["RNF43"].tolist()
# print(lured_rnf43_muts)

count_659 = 0
count_other = 0

for i in lured_rnf43_muts:
    if "G659Vfs*41" in i:
        # print(i)
        count_659 += 1
    else:
        # print("!!!! ", i)
        count_other += 1

print("\n")
print("Number of 659fs:", count_659)
print("Number of other mutations:", count_other)


#no BRAF
lured_rnf43_muts_noBRAF = RNF43_muts[RNF43_muts.index.isin(RNF43_catches_noBRAF)]["RNF43"].tolist()
# print(lured_rnf43_muts)
count_659_nobraf = 0
count_other_nobraf = 0

for i in lured_rnf43_muts_noBRAF:
    if "G659Vfs*41" in i:
        # print(i)
        count_659_nobraf += 1
    else:
        # print("!!!! ", i)
        count_other_nobraf += 1

print("\n")
print("No BRAF- Number of 659fs:", count_659_nobraf)
print("No BRAF- Number of other mutations:", count_other_nobraf)

#saving clinical info for RNF43 catches
clinical[clinical.index.isin(RNF43_catches)].to_csv("outputs/coadread_tcga_clinical_data_RNF43_lured.csv")
clinical[clinical.index.isin(RNF43_catches_noBRAF)].to_csv("outputs/coadread_tcga_clinical_data_RNF43_lured_nobraf.csv")
clinical2[clinical2.index.isin(RNF43_catches)].to_csv("outputs/coadread_tcga_pan_can_atlas_2018_clinical_data_RNF43_lured.csv")
clinical2[clinical2.index.isin(RNF43_catches_noBRAF)].to_csv("outputs/coadread_tcga_pan_can_atlas_2018_clinical_data_RNF43_lured_nobraf.csv")


#######################
# lollipop
#######################

#flatten list of all RNF43 mutations in LURED samples
lured_rnf43_muts_ = [ i.split() for i in lured_rnf43_muts]
lured_rnf43_muts_flat =  [item for sublist in lured_rnf43_muts_ for item in sublist]
count  = Counter(lured_rnf43_muts_flat)

cmd = "../Fig1A-RNF43-lollipops/lollipops-v1.5.2-mac64/lollipops  -domain-labels=off  -w=600 -o outputs/lolli_out.svg RNF43 "

#plot only * truncating
for key, val in count.items():
    if "*"  in key:
        if val == 1:
            size=1
        else:
            size=val**3.5
        if "*" in key:
            color= "#000000"
        else:
            color = "#000000"

        cmd += str(key)+color+"@"+str(size)+" "

print(cmd)
os.system(cmd)
