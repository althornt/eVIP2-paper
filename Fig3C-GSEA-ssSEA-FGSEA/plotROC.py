import pandas as pd
import seaborn as sns
import glob
import matplotlib.pyplot as plt
plt.switch_backend('agg')



def main():

    #import GSEA results
    GSEA_sens = pd.read_csv("GSEA-ssGSEA/significant_combined_GSEA_results_sensitivity.csv", index_col=0)
    GSEA_sens.loc["TPR"] = (GSEA_sens.loc["Sensitivity"])/100
    GSEA_sens.loc["FPR"] = (100- GSEA_sens.loc["Specificity"])/100
    print(GSEA_sens.shape)

    #import FGSEA results
    FGSEA_sens = pd.read_csv("FGSEA/significant_combined_FGSEA_results_sensitivity.csv", index_col=0)
    print(FGSEA_sens)


    plt.scatter(x=1-.875, y=1, color='orange', s = 100, marker = "D")
    sns.scatterplot(data=GSEA_sens.T, x="FPR", y="TPR", alpha=0.5, s = 140, color="grey",marker="o")
    sns.scatterplot(data=FGSEA_sens.T, x="FPR", y="TPR", alpha=0.5, s = 140, color="blue", marker="x")
    # plt.title("Sensitivity and specifiicity across 80 GSEA runs")

    plt.legend(bbox_to_anchor=(1, 1.25 ),borderaxespad=0, labels=["eVIP2",'GSEA',"FGSEA"], fontsize = 12)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")


    plt.savefig("ROC.png", dpi = 400, bbox_inches='tight')
    plt.savefig("ROC.svg", dpi = 400, bbox_inches='tight')
    plt.clf()
    plt.close()







if __name__ == "__main__": main()
