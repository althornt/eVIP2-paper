import pandas as pd
import seaborn as sns
import glob
import matplotlib.pyplot as plt
plt.switch_backend('agg')



def main():

    all_files = glob.glob("/Users/alexis/Desktop/eVIP2-appeal-version/analysis/FGSEA/FGSEA_out/*.csv")
    print(all_files)

    li = []
    for filename in all_files:
        df = pd.read_csv(filename, index_col="pathway")[["padj"]]
        df = df.rename({"padj": filename.split("/")[-1].replace(".csv","")}, axis='columns')
        li.append(df)

    frame = pd.concat(li, axis=1, ignore_index=False,sort=True)
    print(frame)

    frame.to_csv("combined_FGSEA_results.csv")

    sig_frame = frame[frame<.05]

    #calculate sensitivity and specificity
    calc(sig_frame)


def calc(df):
    #subset to the pathways we have validation on
    all_assay_pathways= ["HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_NOTCH_SIGNALING",
    "HALLMARK_P53_PATHWAY","HALLMARK_TGF_BETA_SIGNALING","HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2","HALLMARK_HYPOXIA","HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_KRAS_SIGNALING_DN","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
    positive_assay = ["HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_HYPOXIA",
        "HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
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
    df.loc["Sensitivity"] = df.loc["TP"]/(df.loc["TP"]+df.loc["FN"])
    df.loc["Specificity"] = df.loc["TN"]/(df.loc["FP"]+df.loc["TN"])
    df.loc["PPV"] = df.loc["TP"]/(df.loc["TP"]+df.loc["FP"])
    df.loc["NPV"] = df.loc["TN"]/(df.loc["FN"]+df.loc["TN"])
    df.loc["TPR"] = df.loc["Sensitivity"]
    df.loc["FPR"] = 1- df.loc["Specificity"]

    print(df)
    df.to_csv("significant_combined_FGSEA_results_sensitivity.csv")

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


    #import GSEA results
    GSEA_sens = pd.read_csv("significant_combined_GSEA_results_sensitivity.csv", index_col=0)
    print(GSEA_sens)


    GSEA_sens.loc["TPR"] = (GSEA_sens.loc["Sensitivity"])/100
    GSEA_sens.loc["FPR"] = (100- GSEA_sens.loc["Specificity"])/100

    print(GSEA_sens)


    plt.scatter(x=1-.875, y=1, color='orange', s = 100, marker = "D")
    sns.scatterplot(data=GSEA_sens.T, x="FPR", y="TPR", alpha=0.5, s = 140, color="grey",marker="o")
    sns.scatterplot(data=df.T, x="FPR", y="TPR", alpha=0.5, s = 140, color="blue", marker="x")
    # plt.title("Sensitivity and specifiicity across 80 GSEA runs")

    plt.legend(loc='lower right', labels=["eVIP2",'GSEA',"FGSEA", ], fontsize = 12)

    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")


    plt.savefig("ROC.png", dpi = 400)

    plt.savefig("ROC.svg", dpi = 400)
    plt.clf()
    plt.close()







if __name__ == "__main__": main()
