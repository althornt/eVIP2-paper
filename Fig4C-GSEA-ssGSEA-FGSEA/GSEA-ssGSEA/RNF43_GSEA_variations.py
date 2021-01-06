import pandas as pd
import gseapy
import glob
import sys
from scipy.stats import zscore
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')

"""
GSEA
Two RNF43-G659fs relevant comparisons:
1. "RNF43_G659fs","GFP"
2. "RNF43_G659fs","RNF43_WT"

Four types of files
1. kallisto genes TPM
2. kallisto genes TPM filtered and log transformed
3. z-score on kallisto genes TPM filtered and log transformed
4. mutation specific

Two permutation types
1. gene_set
2. phenotype

Method used to calculate correlation or ranking
1. "signal_to_noise"
2. "t_test"
3. "ratio_of_classes"
4. "diff_of_classes"
5. "log2_ratio_of_classes"

ssGSEA
Four types of files
1. kallisto genes TPM
2. kallisto genes TPM filtered and log transformed
3. z-score on kallisto genes TPM filtered and log transformed
4. mutation specific

Three types of sample_norm_method
1. "rank"
2. "log"
3. "log_rank"


"""

def main():
    # global out_dir
    # out_dir = "/private/groups/brookslab/althornt/eVIP-analysis/eVIP2-revision-analysis-2020/RNF43_GSEA_r_12172020"
    eVIP2_out_dir = "/Users/alexis/Desktop/eVIP2/tutorial_files/eVIP2_out_20201223/"
    kallisto_genes_file = eVIP2_out_dir+"kallisto_files/combined_kallisto_abundance_genes.tsv"

    #for log transformed data, need to replace negative values with 0 for it to run
    filtered = pd.read_csv(eVIP2_out_dir+"kallisto_files/combined_kallisto_abundance_genes_filtered_transformed.tsv", sep="\t", index_col="#gene_id")
    filtered[filtered < 0] = 0
    filtered.to_csv("combined_kallisto_abundance_genes_filtered_transformed_pos.tsv", sep="\t")

    kallisto_genes_filtered_transformed = "combined_kallisto_abundance_genes_filtered_transformed_pos.tsv"

    #calc zscore row-wise
    z_array = zscore(filtered, axis=1, ddof=1)
    z_df = pd.DataFrame(data=z_array, index=filtered.index, columns=filtered.columns)
    z_df.to_csv("combined_kallisto_abundance_genes_filtered_transformed_pos_zscores.tsv", sep="\t")

    zscore_file = "combined_kallisto_abundance_genes_filtered_transformed_pos_zscores.tsv"
    RNF43_G659fs_mutspec = eVIP2_out_dir+"kallisto_files/combined_kallisto_abundance_genes_filtered_transformed_RNF43_G659fs_mutspec.tsv"

    ##################
    # running GSEA
    ##################
    inputs = [kallisto_genes_file, kallisto_genes_filtered_transformed, zscore_file, RNF43_G659fs_mutspec]
    input_label = ["TPM","TPM_filtered_transformed", "TPM_filtered_transformed_zscore", "RNF43_G659fs_mutspec_TPM_filtered_transformed"]
    # comparisons = [["RNF43_G659fs","GFP"],["RNF43_G659fs","RNF43_WT"],["RNF43_R117fs","GFP"],["RNF43_R117fs","RNF43_WT"],["RNF43_WT","GFP"],["RNF43_R117fs","RNF43_G659fs"],["RNF43_G659fs","RNF43_R117fs"]]
    comparisons = [["RNF43_G659fs","GFP"],["RNF43_G659fs","RNF43_WT"]]


    for comparison in comparisons:
        name = ("_vs_").join(comparison)
        print("\n")
        print(name)
        for input,label in zip(inputs,input_label):
            for permtype in ["gene_set","phenotype"]:
                for method in ["signal_to_noise","t_test","ratio_of_classes","diff_of_classes","log2_ratio_of_classes"]:
                    variation = name+"_"+label+"_"+permtype+"_"+method
                    GSEA(comparison,input,variation,permtype,method)


    # comparing all results
    all_files = glob.glob("GSEA_out/*/*.csv")

    li = []
    for filename in all_files:
        df = pd.read_csv(filename, index_col="Term")[["fdr"]]
        df = df.rename({"fdr": filename.split("/")[-2]}, axis='columns')
        li.append(df)

    frame = pd.concat(li, axis=1, ignore_index=False,sort=True)
    frame.to_csv("combined_GSEA_results.csv")

    sig_frame = frame[frame<.25]

    # print("GSEA variations with 0 signifant pathways:")
    # null_variations = sig_frame.columns[sig_frame.isnull().all(0)].tolist()
    # print(len(null_variations))
    # for i in null_variations:
    #     print(i)

    sig_frame = sig_frame.dropna(how='all')
    print(sig_frame)
    sig_frame.to_csv("significant_combined_GSEA_results.csv")

    ##########################
    # Evaluating results
    ###########################
    sig_frame = pd.read_csv("significant_combined_GSEA_results.csv",index_col=0)
    calc(sig_frame)

    ##########################
    # Running SSEA
    ###########################

    ssGSEA_dfs=  []

    for input,label in zip(inputs,input_label):
        print(label)

        df = pd.read_csv(input, sep="\t", index_col="#gene_id")
        filter_col = [col for col in df if col.startswith("RNF43_G659fs")]
        df = df[filter_col]

        for method in  ["rank", "log", "log_rank"]:

            ssGSEA_out = gseapy.ssgsea(data=df, \
                        gene_sets="../h.all.v6.0.symbols.gmt", \
                        outdir="ssGSEA_out/"+label+"_"+method,\
                        processes=20, min_size = 10, \
                        sample_norm_method = method, no_plot=True,seed=7)

            all_assay_pathways= ["HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_NOTCH_SIGNALING",
            "HALLMARK_P53_PATHWAY","HALLMARK_TGF_BETA_SIGNALING","HALLMARK_E2F_TARGETS",
            "HALLMARK_G2M_CHECKPOINT","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_MYC_TARGETS_V1",
            "HALLMARK_MYC_TARGETS_V2","HALLMARK_HYPOXIA","HALLMARK_KRAS_SIGNALING_UP",
            "HALLMARK_KRAS_SIGNALING_DN","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]

            ssGSEA_out = ssGSEA_out.res2d[ssGSEA_out.res2d.index.isin(all_assay_pathways)]
            ssGSEA_dfs.append(ssGSEA_out)

    all_ssGSEA = pd.concat(ssGSEA_dfs,axis=1)
    all_ssGSEA = all_ssGSEA.sort_index()
    print(all_ssGSEA)

    all_ssGSEA.to_csv("ssGSEA_NES.csv")

    sns.stripplot(data=all_ssGSEA.T,edgecolor="black",alpha=.05,size=8,linewidth=1.0,jitter=.01,color="blue")
    plt.title("NES across ssGSEA runs for each RNF43 G659fs sample")
    plt.ylabel("Normalized Enrichment Score")
    plt.xlabel("Pathway")

    plt.xticks(rotation=90, ha='right')
    plt.savefig("ssGSEA_NES_dist.png",bbox_inches = "tight", dpi = 400)
    plt.clf()
    plt.close()


######################################################################
# Functions
######################################################################

def calc(df):
    #subset to the pathways we have validation on
    all_assay_pathways= ["HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_NOTCH_SIGNALING",
    "HALLMARK_P53_PATHWAY","HALLMARK_TGF_BETA_SIGNALING","HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2","HALLMARK_HYPOXIA","HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_KRAS_SIGNALING_DN","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
    positive_assay = ["HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_HYPOXIA","HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
    negative_assay = [i for i in all_assay_pathways if i not in positive_assay]

    df=df[df.index.isin(all_assay_pathways)]

    #need to make a copy to not affect counting of NANs vs vals in df
    df_out = df.copy()

    #TP = positive assay and in significant GSEA df
    df_out.loc["TP"]= df[df.index.isin(positive_assay)].count()
    #TN = not  in positive assay and not significant GSEA
    df_out.loc["TN"] = df[df.index.isin(negative_assay)].isna().sum()
    #FP = negative in assay , significant in GSEA
    df_out.loc["FP"] = df[df.index.isin(negative_assay)].count()
    #FN = positive in assay, and not significant in GSEA
    df_out.loc["FN"] = df[df.index.isin(positive_assay)].isna().sum()

    df = df_out.copy()
    df.loc["Sensitivity"] = df.loc["TP"]/(df.loc["TP"]+df.loc["FN"])*100
    df.loc["Specificity"] = df.loc["TN"]/(df.loc["FP"]+df.loc["TN"])*100
    df.loc["PPV"] = df.loc["TP"]/(df.loc["TP"]+df.loc["FP"])*100
    df.loc["NPV"] = df.loc["TN"]/(df.loc["FN"]+df.loc["TN"])*100

    df.loc["TPR"] = (df.loc["Sensitivity"])/100
    df.loc["FPR"] = (100- df.loc["Specificity"])/100

    df.to_csv("significant_combined_GSEA_results_sensitivity.csv")

    # print(df.T[["Sensitivity","Specificity"]])

    KRASdf  = df[df.index.isin(["HALLMARK_KRAS_SIGNALING_DN","HALLMARK_KRAS_SIGNALING_UP"])]

    """
    #GSEA runs  that dont  call either KRAS pathway
    KRASdf_null = sorted(KRASdf.columns[KRASdf.isnull().all(0)].tolist())
    print(KRASdf_null)
    print(len(KRASdf_null))

    for i in KRASdf_null:
        print(i)

    print("\n")
    #GSEA runs  that  call a KRAS pathway
    KRASdf_notull = sorted(KRASdf.dropna(how="all",axis="columns").columns.tolist())

    print(KRASdf_notull)
    for i in KRASdf_notull:
        print(i)
    """

    #make scatter plot for sensitivity vs specifiicity , which run is the best?
    sns.scatterplot(data=df.T, x="Sensitivity", y="Specificity", alpha=0.5, s = 140)
    plt.scatter(x=100, y=87.5, color='r', s = 100)
    plt.title("Sensitivity and specifiicity across 80 GSEA runs")
    plt.legend(loc='lower left', labels=['GSEA', 'eVIP2'], fontsize = 12)

    plt.savefig("sensitivity_vs_specificity.png", dpi = 400)
    plt.clf()
    plt.close()

def GSEA(samples,file_path,out_dir_name,permtype,method):
    """
    samples: the names of the 2 conditions to compare
    file_path: location of gene expression input
    out_dir_name: what to name the output directory
    permtype: type of permutation
    """

    print("----------------------------------------------------------------------------------------")
    print(out_dir_name)
    print("\n")

    hallmark_gmt = "../h.all.v6.0.symbols.gmt"

    #setting conditions
    samples_cls = [samples[0]]*4 + [samples[1]]*4

    #subset df to give samples
    df = pd.read_csv(file_path, sep="\t", index_col="#gene_id")
    filter_col = [col for col in df if col.startswith(tuple(samples))]
    df_filter=df[filter_col]
    # print(df_filter.head())

    #running GSEA
    gseapy.gsea(data=df_filter, \
                gene_sets=hallmark_gmt, \
                cls=samples_cls, outdir="GSEA_out/"+out_dir_name,\
                processes=20, permutation_type = permtype, \
                min_size = 10, method = method, no_plot=True, seed=123)


    #signficant results (FDR < .25)
    report = pd.read_csv("GSEA_out/"+out_dir_name+"/gseapy.gsea."+permtype+".report.csv",index_col='Term')
    print(report[report['fdr']<.25])


if __name__ == "__main__": main()
