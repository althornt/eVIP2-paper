import pandas as pd
from collections import Counter
import os

#files from cbioportal RNF43 mutation downloads

def format(file):
    print("\n")
    print(file)
    df = pd.read_csv(file, sep="\t")

    l = df["RNF43"].tolist()
    l = [i.split() for i in l]
    flat_list = [item for sublist in l for item in sublist if item not in ["WT","NS"] ]
    li = [i.split("*")[0] for i in flat_list ]

    #remove  extra character before fs so in order to combine fs the same position
    new_list = []
    for  i in li:
        if i.endswith("fs"):
            new_list.append(i[:4]+i[5:])
        else:
            new_list.append(i)

    count  = Counter(new_list)
    print(count)
    cmd = "./lollipops-v1.5.2-mac64/lollipops -domain-labels=off -legend -show-motifs -w=1150 -o "+file+"_out.svg  RNF43 "

    for key, val in count.items():
        # print(key,val)

        if val == 1:
            size=.01
        else:
            size=val**3.5

        if key.split("*")[0].endswith("fs"):
            color= "#228b22"
        else:
            color = "#000000"

        cmd += str(key)+color+"@"+str(size)+" "

    print(cmd)

    return(cmd)


CRC_DFCI_TCGA_comd = format("RNF43_mutations_Colorectal_adenocarcinoma_DFCI_and_TCGApancancer.txt")
os.system(CRC_DFCI_TCGA_comd)

Endo_TCGA_comd = format("RNF43_mutations_Uterine_Corpus_Endometrial_Carcinoma_TCGAPanCancer.txt")
os.system(Endo_TCGA_comd)

Stom_TCGA_comd = format("RNF43_mutations_Stomach_Adenocarcinoma_TCGAPanCancer.txt")
os.system(Stom_TCGA_comd)


# CRC_TCGA_comd = format("ColorectalCancerTCGA.txt")
# os.system(CRC_TCGA_comd)
#
# Endo_TCGA_comd = format("EndometrialCancerTCGA.txt")
# os.system(Endo_TCGA_comd)
#
# CRC_NHS_comd = format("ColorectalCancerNHS-HPFS.txt")
# os.system(CRC_NHS_comd)



    # df_counts = df["Protein_Change"].value_counts()
    #
    #
    # cmd = "./lollipops-v1.5.2-mac64/lollipops -legend -o "+file+"_out.svg RNF43 "
    # # cmd += file +"_out"+" "
    #
    # print(df_counts)
    #
    # for index, val in df_counts.iteritems():
    #     if val == 1:
    #         size=.1
    #     else:
    #         size=val**2
    #
    #     if index.endswith("fs"):
    #         color= "#228b22"
    #     else:
    #         color = "#000000"
    #
    #     # if index in ["p.G659fs","p.R117fs"]:
    #     #     cmd += str(index)+"#00ff00"+"@"+str(val)+" "
    #     #
    #     # else:
    #     #     cmd += str(index)+"#000000"+"@"+str(val)+" "
    #
    #     cmd += str(index)+color+"@"+str(size)+" "
    #
    #
    # print(cmd)
    # return(cmd)
