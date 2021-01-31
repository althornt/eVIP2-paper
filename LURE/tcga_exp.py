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
import json

def main():
    #from LURE_lolli_
    score_mut = pd.read_csv("outputs/LUREscore_mut_status.csv", index_col=0)

    #Xena - gene-leveltranscription estimates, as in log2(x+1) transformed RSEM normalized count
    # tcga_exp = pd.read_csv("inputs/HiSeqV2", sep="\t",index_col=0)
    tcga_exp_coadread = pd.read_csv("inputs/XenaCOADREAD_HiSeqV2", sep="\t",index_col=0)
    tcga_exp_coadread.columns =  [i[:-3] for i in tcga_exp_coadread.columns]

    #add  kras status to mut df
    kras_mut =  pd.read_csv("inputs/cbioportal_KRAS_TCGA_COAD_mutations.txt", sep="\t",index_col  = "SAMPLE_ID")
    kras_mut.index =  [i[:-3] for i in kras_mut.index]

    score_mut = score_mut.merge(kras_mut["KRAS"], how='outer', left_index=True, right_index=True)
    score_mut = score_mut[score_mut['x'].notna()]

    ##########################
    # get gene lists
    #########################

    # #import kegg gene set
    # kegg_mapk_genes = None
    # wnt_genes = None
    # cancer_genes = None
    # cytokine_genes = None
    # all_kegg_genes = []
    # with open("inputs/c2.cp.kegg.v7.2.symbols.gmt", "r") as kegg:
    #     for line in kegg:
    #         all_kegg_genes.extend(line.split()[2:])
    #         if line.split()[0] =="KEGG_MAPK_SIGNALING_PATHWAY":
    #             kegg_mapk_genes = line.split()[2:]
    #         if line.split()[0] =="KEGG_WNT_SIGNALING_PATHWAY":
    #             wnt_genes = line.split()[2:]
    #         if line.split()[0] =="KEGG_PATHWAYS_IN_CANCER":
    #             cancer_genes = line.split()[2:]
    #         if line.split()[0] =="KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION":
    #             cytokine_genes = line.split()[2:]


    # nfkb_genes = []
    # with open("inputs/h.all.v6.0.symbols.gmt", "r") as hm:
    #     for line in hm:
    #         # if line.split()[0] in ["HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_HYPOXIA","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]:
    #         if line.split()[0] in ["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]:
    #
    #         # if line.split()[0] in ["HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN"]:
    #             nfkb_genes.extend(line.split()[2:])


    #https://www.nature.com/articles/s41698-018-0051-4 "MAPK Pathway Activity Score "
    global mapkgenes
    mapkgenes = "PHLDA1, SPRY2, SPRY4, ETV4, ETV5, DUSP4, DUSP6, CCND1, EPHA2,  EPHA4".replace(',','').split()
    # wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3"]
    # wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3","DKK1","WIF1","LRP6","CDH1","CTNNB1","GSK3B"]
    global wnt
    wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3","WIF1"]

    # wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3","DKK1","WIF1"]

    global  mmr
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

    # clustermap(score_mut,tcga_exp_coadread,kegg_mapk_genes+mapkgenes+wnt_genes+mmr,"kegg_mapk_genes_andtengene" )

    lobodo_ras = open("inputs/lobodo_RAS_signature.txt","r").read().splitlines()
    # clustermap(score_mut,tcga_exp_coadread,lobodo_ras,"lobodo_ras" )

    popovici = pd.read_csv( "inputs/popovici_braf_signature.tsv",sep="\t").drop(columns=["Pair","Pair.1"])
    popovici_genes = [item for sublist in popovici.values.tolist() for item in sublist]

    # clustermap(score_mut,tcga_exp_coadread,mapkgenes,"mapk" )
    # clustermap(score_mut,tcga_exp_coadread,popovici_genes+mapkgenes,"braf_sig_mapk" )
    # clustermap(score_mut,tcga_exp_coadread,mapkgenes+wnt+mmr,"og" )


    ##########################
    # make clustermaps
    #########################


    with open("inputs/hallmark.all.v6.0_c2.cp.kegg.v6.0_c2.cp.reactome.v6.0_symbols.json") as f:
        pathway_dict = json.load(f)
    all_pathway_genes =  [item for sublist in  pathway_dict.values() for item in sublist]

    """
    dusp4_like_df = tcga_exp_coadread.loc[set(all_pathway_genes)].transpose().corr()["DUSP4"].sort_values(ascending=False)
    dusp4_like_head = dusp4_like_df.head(30).index.tolist()
    dusp4_like_tail = dusp4_like_df.dropna().tail(30).index.tolist()
    print(dusp4_like_head)
    print(dusp4_like_tail)

    etv5_like_df = tcga_exp_coadread.loc[set(all_pathway_genes)].transpose().corr()["ETV5"].sort_values(ascending=False)
    etv5_like_head = etv5_like_df.head(30).index.tolist()
    etv5_like_tail = etv5_like_df.dropna().tail(30).index.tolist()
    print(etv5_like_head)
    print(etv5_like_tail)
    """

    pathways_eVIP =["HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            "HALLMARK_HYPOXIA","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","KEGG_MAPK_SIGNALING_PATHWAY",
            "REACTOME_MAPK_TARGETS_NUCLEAR_EVENTS_MEDIATED_BY_MAP_KINASES","REACTOME_ERK_MAPK_TARGETS"]

    eVIP_pway_genes = []
    for i in pathways_eVIP:
        eVIP_pway_genes.extend(pathway_dict[i])

    # for i in pathways_eVIP:
    #     for j in pathway_dict[i]:
    #         if j in dusp4_like:
    #             print(i,j)
    #         if j in etv5_like:
    #             print(i,j)

    # pathway_matches = {}
    #
    # for i in [dusp4_like+etv5_like]:
    #     for k,v in pathway_dict.items():
    #         if i in v:
    #             print(i,k)
    #             try:
    #                 pathway_matches[k].append(i)
    #             except:
    #                 pathway_matches = [i]

    # from collections import defaultdict
    #
    # pathway_matches_dict = defaultdict(list)
    # for i in dusp4_like+etv5_like:
    #     for k,v in pathway_dict.items():
    #         if i in v:
    #             print(i)
    #             pathway_matches_dict[k].append(i)
    #
    #
    # print(pathway_matches_dict)
    # for k in sorted(pathway_matches_dict, key=lambda k: len(pathway_matches_dict[k]), reverse=True):
    #     if k.startswith("RE"):
    #         print(k,v)
    #

    """
    ############################################################
    # Finding high corr genes in RNF43 trunc samples
    ############################################################

    #find tcga samples in both df
    samples = [i for i in tcga_exp_coadread.columns if i in score_mut.index ]
    tcga_exp = tcga_exp_coadread[samples]
    score_mut_rnf43 = score_mut[score_mut.index.isin(samples)]["RNF43_TRUNCATING"].to_frame()
    tcga_exp = tcga_exp.append(score_mut_rnf43["RNF43_TRUNCATING"].transpose())

    #get samples with RNF43 trunc
    tcga_exp_rnf43 = tcga_exp.loc[:,tcga_exp.loc["RNF43_TRUNCATING"]=="TRUNCATING"].drop(["RNF43_TRUNCATING"], axis=0)

    #reduce to eVIP pway genes
    tcga_exp_rnf43 = tcga_exp_rnf43[tcga_exp_rnf43.index.isin(eVIP_pway_genes)]

    #find correlated genes
    corr = tcga_exp_rnf43.astype(float).transpose().corr().abs()
    corr = corr.unstack()
    corr = corr.sort_values(kind="quicksort",ascending=False)
    corr_clean = corr.dropna()
    corr_clean = corr_clean[corr_clean>0.95]
    corr_clean = corr_clean[corr_clean<0.999999999] #remove self-correlations
    candidate_genes = corr_clean.index.tolist()
    candidate_genes = list(set([item for sublist in candidate_genes for item in sublist]))
    print(len(candidate_genes))
    """
    """
    hypoxia tried - VHL     HIF1AN

    EMT - CDH1 is down in half but makes clustering worse
        - VIM is up in  half but makes clustering worse (same half as CDH1)
        FN1
        TJP1
        OCLN
        CDH2
        DKK1
    """



    global hypoxia_genes,emt_genes,nfkb_genes
    hypoxia_genes = ["HIF1A"]

    # global
    emt_genes = ["TWIST1","SNAI1"]

    # global nfkb_genes
    nfkb_genes = ["NFKB1","NFKB2","RELB"]


    for df in [tcga_exp_coadread]:
        clustermap(score_mut,df,mapkgenes+wnt+mmr+hypoxia_genes+nfkb_genes+emt_genes,"og" )





def clustermap(score_mut,tcga_exp, gene_list,  label ):
    tcga_exp = tcga_exp[tcga_exp.index.isin(gene_list)]

    #remove low std genes
    # tcga_exp = tcga_exp.loc[tcga_exp.std(axis=1) > 2, :]

    #drop duplicated samples
    tcga_exp = tcga_exp.loc[:,~tcga_exp.columns.duplicated(keep="first")]

    #find tcga samples in both
    samples = [i for i in tcga_exp.columns if i in score_mut.index ]
    tcga_exp = tcga_exp[samples]
    score_mut = score_mut[score_mut.index.isin(samples)][["RNF43_TRUNCATING","BRAF_MISSENSE","x","RNF43","KRAS"]]
    #rename RNF43 so its not confused with gene expression row
    score_mut = score_mut.rename(columns={"RNF43": "RNF43_muts","KRAS":"KRAS_muts"})
    tcga_exp = tcga_exp.append(score_mut[["RNF43_TRUNCATING","BRAF_MISSENSE","x","RNF43_muts","KRAS_muts"]].transpose())

    xenaclinical = pd.read_csv("inputs/COADREAD_clinicalMatrix", sep="\t",
                            index_col = "_PATIENT")
    xenaclinical = xenaclinical[xenaclinical.index.isin(tcga_exp.columns)]
    xenaclinical = xenaclinical[xenaclinical.index.isin(tcga_exp)]
    xenaclinical = xenaclinical.loc[~xenaclinical.index.duplicated(keep="first")]
    tcga_exp = tcga_exp.append(xenaclinical["CDE_ID_3226963"].transpose())

    #adding driver status of all mapk/rtk genes
    mapk_driver_status = pd.read_csv("inputs/cbioportal-MAPK-RTK-genes-driver-sample-matrix.txt",
                    sep="\t", index_col = 0)
    mapk_driver_status.index =  [i.split(":")[1][:-3] for i in mapk_driver_status.index]
    mapk_driver_status = mapk_driver_status[mapk_driver_status.index.isin(tcga_exp.columns)]

    tcga_exp = tcga_exp.append(mapk_driver_status["Altered"].transpose())

    ##################################
    # adding color columns
    ##################################

    rnf43_status = tcga_exp.loc["RNF43_TRUNCATING"].tolist()
    braf_status = tcga_exp.loc["BRAF_MISSENSE"].tolist()
    lure_score = tcga_exp.loc["x"].tolist()
    RNF43_muts  = tcga_exp.loc["RNF43_muts"].tolist()
    KRAS_muts = tcga_exp.loc["KRAS_muts"].tolist()
    msi_status = tcga_exp.loc["CDE_ID_3226963"].tolist()
    mapk_pway_driver_status = tcga_exp.loc["Altered"].tolist()

    #convert MAPK driver status to color
    mapk_pway_driver_status_c = []
    for i in mapk_pway_driver_status:
        if i == 0:
            mapk_pway_driver_status_c.append("#f0edf2")
        elif i == 1:
            mapk_pway_driver_status_c.append("black")
        else:
            mapk_pway_driver_status_c.append("#f0edf2")

    #convert MSI status to color
    msi_status_c = []
    for i in msi_status:
        if i == "MSI-H":
            msi_status_c.append("darkred")
        elif i == "MSI-L":
            msi_status_c.append("indianred")
        elif i == "MSS":
            msi_status_c.append("mistyrose")
        else:
            msi_status_c.append("#f0edf2")

    #convert RNF43 mut status to color
    kras_status_c = []
    g12 = ["G12D","G12S", "G12R", "G12C", "G12V", "G12A"]
    for i in KRAS_muts:
        res = any(ele in str(i) for ele in g12)
        if res ==  True:
            kras_status_c.append("gray")
        else:
            kras_status_c.append("#f0edf2")

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


    #prep df for clustermap
    tcga_exp = tcga_exp.drop(["RNF43_TRUNCATING","BRAF_MISSENSE","x","x_color",
                "RNF43_muts","KRAS_muts","CDE_ID_3226963","Altered"], axis=0)
    tcga_exp = tcga_exp.convert_objects(convert_numeric=True)

    #now with only genes rows left, label colors
    gene_colors = []
    for i in tcga_exp.index.tolist():
        if i in mapkgenes:
            gene_colors.append("black")
        elif i in wnt:
            gene_colors.append("cornflowerblue")
        elif i in mmr:
            gene_colors.append("plum")
        elif i in nfkb_genes:
            gene_colors.append("purple")
        elif i in hypoxia_genes:
            gene_colors.append("crimson")
        elif i in emt_genes:
            gene_colors.append("slateblue")
        else:
            gene_colors.append("pink")

    ################
    # clustermap
    ################
    sns.clustermap(tcga_exp,z_score=0,center=0,cmap="bwr",xticklabels=False,
                yticklabels=True,row_colors=gene_colors, method = "ward",
                colors_ratio = .03, figsize = (12,8),
                col_colors=[lure_scale_c,braf_status_c,rnf43_status_c,RNF43_mut_type_c,mapk_pway_driver_status_c,kras_status_c,msi_status_c],
                row_cluster=True)

    plt.tight_layout()
    plt.savefig("outputs/clustermaps/tcga_heatmap_"+label+".png",bbox_inches = "tight",dpi=400)
    plt.savefig("outputs/clustermaps/tcga_heatmap_"+label+".svg",bbox_inches = "tight",dpi=400)

    plt.clf()
    plt.close()

    ##################################################################
    # Putting legends and colorbars in another file
    ##################################################################

    #gene row patches
    mapk_patch = mpatches.Patch(color='black', label='MAPK')
    nfkb_patch = mpatches.Patch(color='purple', label='NFKB')
    hypoxia_patch = mpatches.Patch(color='crimson', label='Hypoxia')
    emt_patch = mpatches.Patch(color='slateblue', label='EMT')
    wnt_patch = mpatches.Patch(color='cornflowerblue', label='WNT')
    mmr_patch = mpatches.Patch(color='plum', label='MMR')
    plt.legend(handles=[mapk_patch,nfkb_patch,hypoxia_patch,wnt_patch,emt_patch, mmr_patch],loc='upper center')

    #msi status patch
    msi_high_patch = mpatches.Patch(color='darkred', label='MSI-H')
    msi_low_patch = mpatches.Patch(color='indianred', label='MSI-L')
    mss_patch = mpatches.Patch(color='mistyrose', label='MSS')
    plt.legend(handles=[msi_high_patch,msi_low_patch,mss_patch],loc='lower left')

    #color baar for LURE classifier score
    sm = plt.cm.ScalarMappable(cmap=plt.cm.YlGn, norm=plt.Normalize(vmin=0, vmax=1))
    sm._A = []
    plt.colorbar(sm)

    #color baar for heatmap zscore vaalue
    sm = plt.cm.ScalarMappable(cmap=plt.cm.bwr, norm=plt.Normalize(vmin=-5, vmax=5))
    sm._A = []
    plt.colorbar(sm)

    plt.savefig("outputs/tcga_heatmap_legend.png",bbox_inches = "tight",dpi=400)
    plt.savefig("outputs/tcga_heatmap_legend.svg",bbox_inches = "tight",dpi=400)
    plt.clf()
    plt.close()


if __name__ == "__main__":
    main()
