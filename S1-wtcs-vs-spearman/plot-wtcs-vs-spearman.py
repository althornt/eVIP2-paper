import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy import stats

"""
Comparing adjusted eVIP p-values when using different correlation measures
(L1000 weighted connectivity score vs. spearman correlation).

eVIP commands:

#using z-score only, eVIP calculates spearman
python /Users/alexis/Desktop/eVIP2/run_eVIP.py
-zscore_gct high_rep_A549_8reps_ZSPCLM_150703_n4232x978.gct
-out_directory /Users/alexis/Desktop/eVIP-L1000/zscore_eVIP_out
-sig_info high_rep_A549_8reps_COMPZ.MODZ_150703_sig.info
-c lung_control_x_mutation_status.grp
-r TA_clone_info_MASTER_20150514_150629_download_WTref_x_mut_fixed_150706_INCLUDE_only_ref2allele.txt
-allele_col x_mutation_status -ie_col x_ie_a549 -num_reps 8 -x_thresh 1.3
-y_thresh 1.3 -xmax 7 -ymin -4.5 -ymax 4 -use_c_pval -i 1000 -ie_filter 0.5
-mut_wt_rep_thresh 0.05 -mut_wt_rep_rank_diff 0 -disting_thresh 0.05


#using z-score and wtcs
python /Users/alexis/Desktop/eVIP2/run_eVIP.py -wtcs_gct self_score_n4232x4232.gct
-zscore_gct high_rep_A549_8reps_ZSPCLM_150703_n4232x978.gct
-out_directory /Users/alexis/Desktop/eVIP-L1000/wtcs_eVIP_out
-sig_info high_rep_A549_8reps_COMPZ.MODZ_150703_sig.info
-c lung_control_x_mutation_status.grp
-r TA_clone_info_MASTER_20150514_150629_download_WTref_x_mut_fixed_150706_INCLUDE_only_ref2allele.txt
-allele_col x_mutation_status -ie_col x_ie_a549 -num_reps 8 -x_thresh 1.3
-y_thresh 1.3 -xmax 7 -ymin -4.5 -ymax 4 -use_c_pval -i 1000 -ie_filter 0.5
-mut_wt_rep_thresh 0.05 -mut_wt_rep_rank_diff 0 -disting_thresh 0.05

"""

#import  predict file using z-score and WTCS
wtcs = pd.read_csv("wtcs_predict.txt",
                    sep="\t", index_col  = "mut")

#import predict file from using z-score spearman
spearman = pd.read_csv("spearman_predict.txt",
                    sep="\t", index_col  = "mut")

#supplement from Berger et al
paper = pd.read_csv("eVIP_paper_L1000_calls.txt",
                    sep = "\t", index_col = "Allele")["eVIP_prediction"]

df_pred =  pd.merge(wtcs["prediction"].rename("wtcs_pred"),
            spearman["prediction"].rename("spearman_pred"),
            left_index=True, right_index=True)

df_pval = pd.merge(wtcs["wt_mut_rep_vs_wt_mut_conn_c_pval"].rename("wtcs_pval"),
            spearman["wt_mut_rep_vs_wt_mut_conn_c_pval"].rename("spearman_pval"),
            left_index=True, right_index=True)

df_pval = -np.log10(df_pval)


pred_pval =  pd.merge(df_pval, df_pred, left_index=True, right_index=True)
print(pred_pval.head())

#verifying wtcs matches Berger et al
pred_pval_paper =  pd.merge(pred_pval, paper, left_index=True, right_index=True)
print((str((pred_pval_paper["wtcs_pred"] ==
    pred_pval_paper["eVIP_prediction"]).sum())) +" of "+
    str(pred_pval_paper.shape[0])+ " match Berger et al")

color_dict={"COF":"#6a3d9a",
"GOF":"#ca0020",
"LOF":"#0571b0",
"Neutral":"black",
"NI":"grey"}

fig = plt.gcf()
fig.set_size_inches(8, 8)

#plotting on top of eachother to create edges,
#since edgecolor cant be used with the pallete dict

#edge color is from spearman
sns.scatterplot(x='wtcs_pval', y='spearman_pval', data=pred_pval, hue = "spearman_pred",
                        palette=color_dict, s  = 50, edgecolor=None)
#middle fill is from WTCS
sns.scatterplot(x='wtcs_pval', y='spearman_pval', data=pred_pval, hue = "wtcs_pred",
                        palette=color_dict,  s = 15 ,edgecolor=None)

#Plot the regression line
sns.regplot(x='wtcs_pval', y='spearman_pval', data=pred_pval, scatter=False, color = "black")

#add lines for significance
plt.axvline(x=-np.log10(.05), color = 'grey', linestyle = '--')
plt.axhline(y=-np.log10(.05), color = 'grey', linestyle = '--')

#calculating pearson and rounding
r2 = str(stats.pearsonr(pred_pval['wtcs_pval'], pred_pval['spearman_pval'])[0])[:5]

#box style
props = dict(boxstyle='square', facecolor='grey', alpha=0.1)
plt.text(1.5, 4, "pearson="+r2, fontsize=12,bbox=props)

plt.title('Adjusted p-values using L1000 weighted connectivity score vs spearman rank correlation')
plt.xlabel('-log10(pval) using weighted connectivity score')
plt.ylabel('-log10(pval) using spearman rank correlation')
plt.legend(loc='center right', bbox_to_anchor=(1.4, 0.5), ncol=1)
plt.savefig("wtcs-vs-spearman.png", dpi = 400, bbox_inches = "tight")
plt.savefig("wtcs-vs-spearman.svg", dpi = 400, bbox_inches = "tight")
plt.clf()
plt.close()
