import os
import seaborn as sns
import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

"""
This script clusters TCGA COAD expression data
"""

def main():
    #BRAF-scores and RNF43 mutants generated from LURE lolli script
    score_mut = pd.read_csv("outputs/LUREscore_mut_status.csv", index_col=0)

    #Xena - gene-level transcription estimates, as in log2(x+1) transformed RSEM normalized count
    # tcga_exp_coadread = pd.read_csv("inputs/XenaCOADREAD_HiSeqV2", sep="\t",index_col=0)

    #cbpioportal RSEM
    tcga_exp_coadread = pd.read_csv("inputs/cbioportal_data_RNA_Seq_v2_expression_median.txt",
                                        sep="\t",index_col=0)
    tcga_exp_coadread =  tcga_exp_coadread.drop(columns="Entrez_Gene_Id")
    tcga_exp_coadread.columns =  [i[:-3] for i in tcga_exp_coadread.columns]
    tcga_exp_coadread = np.log2(tcga_exp_coadread+1) #log2transform

    #add kras status to mut df
    kras_mut =  pd.read_csv("inputs/cbioportal_KRAS_TCGA_COAD_mutations.txt",
                                        sep="\t",index_col  = "SAMPLE_ID")
    kras_mut.index =  [i[:-3] for i in kras_mut.index]
    score_mut = score_mut.merge(kras_mut["KRAS"], how='outer',
                                            left_index=True, right_index=True)
    score_mut = score_mut[score_mut['x'].notna()]


    #z-score expression data
    tcga_exp_coadread_zscore = stats.zscore(tcga_exp_coadread, axis=1)
    tcga_exp_coadread = pd.DataFrame(tcga_exp_coadread_zscore,
                columns = tcga_exp_coadread.columns,
                index = tcga_exp_coadread.index)

    ##########################
    # make clustermap
    #########################

    global mapkgenes,wnt,hypoxia_genes,emt_genes,nfkb_genes,mmr
    hypoxia_genes = ["HIF1A"]
    emt_genes = ["TWIST1","SNAI1"]
    nfkb_genes = ["NFKB1","NFKB2","RELB"]
    #https://www.nature.com/articles/s41698-018-0051-4 "MAPK Pathway Activity Score
    mapkgenes = "PHLDA1, SPRY2, SPRY4, ETV4, ETV5, DUSP4, DUSP6, CCND1,\
                            EPHA2, EPHA4".replace(',','').split()
    wnt = ["AXIN2", "NKD1", "RNF43", "ZNRF3","WIF1"]
    mmr = ["MLH1"]

    clustermap(score_mut,tcga_exp_coadread,mapkgenes+wnt+mmr+hypoxia_genes+
                                                nfkb_genes+emt_genes,"full" )

    ############################################
    #plot only WT BRAF and RNF43 TRUNC samples

    score_mut_BRAFWT_RNF43_trunc = score_mut[(score_mut["BRAF_MISSENSE"].isna())
                                                & (~score_mut["RNF43"].isna())]
    clustermap(score_mut_BRAFWT_RNF43_trunc,tcga_exp_coadread,mapkgenes+wnt+
                mmr+hypoxia_genes+nfkb_genes+emt_genes,"_BRAFWT_RNF43_TRUNC" )
    clustermap(score_mut_BRAFWT_RNF43_trunc,tcga_exp_coadread,
                mapkgenes,"_BRAFWT_RNF43_TRUNC_mapk" )

    ############################################
    #plot only WT BRAF and only RNF43 G659fs samples (no other RNF43 TRUNC variants either)

    score_mut_BRAFWT_RNF43_659 = score_mut[(score_mut["BRAF_MISSENSE"].isna()) &
                                (score_mut["RNF43"] == "['G659fs']")]
    clustermap(score_mut_BRAFWT_RNF43_659,tcga_exp_coadread,mapkgenes,
                "_BRAFWT_RNF43_G659fs_mapk" )
    clustermap(score_mut_BRAFWT_RNF43_659,tcga_exp_coadread,mapkgenes+wnt+mmr+
                    hypoxia_genes+nfkb_genes+emt_genes,"_BRAFWT_RNF43_G659fs" )


def clustermap(score_mut,tcga_exp, gene_list,  label ):
    tcga_exp = tcga_exp[tcga_exp.index.isin(gene_list)]

    #remove low std genes
    # tcga_exp = tcga_exp.loc[tcga_exp.std(axis=1) > 2, :]

    #drop duplicated samples
    tcga_exp = tcga_exp.loc[:,~tcga_exp.columns.duplicated(keep="first")]

    #find tcga samples in both df
    samples = [i for i in tcga_exp.columns if i in score_mut.index ]
    tcga_exp = tcga_exp[samples]
    score_mut = score_mut[score_mut.index.isin(samples)][["RNF43_TRUNCATING",
                    "BRAF_MISSENSE","x","RNF43","KRAS"]]
    #rename RNF43 and KRAS so its not confused with gene expression row
    score_mut = score_mut.rename(columns={"RNF43": "RNF43_muts","KRAS":"KRAS_muts"})
    tcga_exp = tcga_exp.append(score_mut[["RNF43_TRUNCATING","BRAF_MISSENSE",
                    "x","RNF43_muts","KRAS_muts"]].transpose())

    #import clinical data
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

    #convert RNF43 mut status to color
    RNF43_mut_type_c = []
    for i in RNF43_muts:
        if "G659fs" in str(i):
            RNF43_mut_type_c.append("gold")
        else:
            RNF43_mut_type_c.append("#f0edf2")

    #adding scale for LURE classifier score
    norm = matplotlib.colors.Normalize(vmin=0, vmax=1, clip=True)
    mapper = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.YlGn)
    tcga_exp.loc['x_color'] = tcga_exp.loc['x'].apply(
                            lambda x: mcolors.to_hex(mapper.to_rgba(x)))
    lure_scale_c = tcga_exp.loc['x_color'].tolist()

    """
    #left vs right side
    clinical = pd.read_csv("inputs/coadread_tcga_clinical_data.tsv", sep="\t",
                            index_col = "Patient ID")

    left = ['Splenic Flexure','Descending Colon','Sigmoid Colon',
                        'Rectosigmoid Junction','Rectum']
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
    sns.clustermap(tcga_exp,center=0,cmap="bwr",xticklabels=False,
                yticklabels=True,row_colors=gene_colors, method = "ward",
                colors_ratio = .03, figsize = (12,8),
                col_colors=[lure_scale_c,braf_status_c,rnf43_status_c,
                            RNF43_mut_type_c,mapk_pway_driver_status_c,
                            kras_status_c,msi_status_c],
                row_cluster=True,vmin=-5, vmax=5)

    plt.tight_layout()
    plt.savefig("outputs/clustermaps/tcga_heatmap_"+label+".png",
                bbox_inches = "tight",dpi=400)
    plt.savefig("outputs/clustermaps/tcga_heatmap_"+label+".svg",
                bbox_inches = "tight",dpi=400)

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
    plt.legend(handles=[mapk_patch,nfkb_patch,hypoxia_patch,wnt_patch,emt_patch,
                mmr_patch],loc='upper center')
    plt.savefig("outputs/clustermaps/tcga_heatmap_legend1.svg",
                bbox_inches = "tight",dpi=400)
    plt.savefig("outputs/clustermaps/tcga_heatmap_legend1.png"
                ,bbox_inches = "tight",dpi=400)
    plt.clf()
    plt.close()

    #msi status patch
    msi_high_patch = mpatches.Patch(color='darkred', label='MSI-H')
    msi_low_patch = mpatches.Patch(color='indianred', label='MSI-L')
    mss_patch = mpatches.Patch(color='mistyrose', label='MSS')
    plt.legend(handles=[msi_high_patch,msi_low_patch,mss_patch],loc='lower left')

    #color bar for LURE classifier score
    sm = plt.cm.ScalarMappable(cmap=plt.cm.YlGn, norm=plt.Normalize(vmin=0, vmax=1))
    sm._A = []
    plt.colorbar(sm)

    #color baar for heatmap zscore vaalue
    sm = plt.cm.ScalarMappable(cmap=plt.cm.bwr, norm=plt.Normalize(vmin=-5, vmax=5))
    sm._A = []
    plt.colorbar(sm)

    plt.savefig("outputs/clustermaps/tcga_heatmap_legend2.png",bbox_inches = "tight",dpi=400)
    plt.savefig("outputs/clustermaps/tcga_heatmap_legend2.svg",bbox_inches = "tight",dpi=400)
    plt.clf()
    plt.close()



if __name__ == "__main__":
    main()
