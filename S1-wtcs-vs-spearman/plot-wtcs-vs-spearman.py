import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from scipy import stats


"""
Comparing adjusted eVIP p-values when using different correlation measures (L1000 weighted connectivity score vs. spearman correlation).

eVIP commands:

### Running eVIP using weighted connectivity score file
*run_eVIP does not support import of correlation matrix so have to run scripts individually

#### eVIP compare
`python /private/groups/brookslab/althornt/apps/eVIP/bin/eVIP_compare.py --sig_info /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/high_rep_A549_8reps_COMPZ.MODZ_150703_sig.info --gctx /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/self_score_n4232x4232.gct -o /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/self_score_eVIP/eVIP_compare_out -r /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/TA_clone_info_MASTER_20150514_150629_download_WTref_x_mut_fixed_150706_INCLUDE_only_ref2allele.txt -c /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/lung_control_x_mutation_status.grp -i 1000 --conn_null /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/high_rep_A549_8reps_150703_eVIP_compare_0.5_IE_thresh_closed_only_150706_conn_null.txt --ie_col x_ie_a549 --ie_filter 0.5 --num_reps 8`

#### eVIP predict
`python /private/groups/brookslab/althornt/apps/eVIP/bin/eVIP_predict.py -i /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/self_score_eVIP/eVIP_compare_out.txt -o /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/self_score_eVIP/eVIP_predict_out --mut_wt_rep_thresh 0.05 --mut_wt_rep_rank_diff 0 --disting_thresh 0.05 --use_c_pval`

##############################

#Running eVIP using z-score file

`python /private/groups/brookslab/althornt/apps/eVIP/run_eVIP.py -zscore_gct /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/high_rep_A549_8reps_ZSPCLM_150703_n4232x978.gct -out_directory /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/zscore_eVIP -sig_info /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/high_rep_A549_8reps_COMPZ.MODZ_150703_sig.info -c /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/lung_control_x_mutation_status.grp -r /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/TA_clone_info_MASTER_20150514_150629_download_WTref_x_mut_fixed_150706_INCLUDE_only_ref2allele.txt -num_reps 8 -ie_filter 0.5  -ie_col x_ie_a549 -i 1000 -conn_null /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/high_rep_A549_8reps_150703_eVIP_compare_0.5_IE_thresh_closed_only_150706_conn_null.txt -use_c_pval -cell_id A549 -plate_id TA.OE014.OE015`

* got corr from run eVIP

#compare
`python /private/groups/brookslab/althornt/apps/eVIP/bin/eVIP_compare.py --sig_info /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/high_rep_A549_8reps_COMPZ.MODZ_150703_sig.info --gctx /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/zscore_eVIP/eVIP_out/spearman_rank_matrix.gct -o /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/zscore_eVIP/eVIP_compare_out -r /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/TA_clone_info_MASTER_20150514_150629_download_WTref_x_mut_fixed_150706_INCLUDE_only_ref2allele.txt -c /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/lung_control_x_mutation_status.grp -i 1000 --conn_null /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/FINAL_RUN_FILES/high_rep_A549_8reps_150703_eVIP_compare_0.5_IE_thresh_closed_only_150706_conn_null.txt --ie_col x_ie_a549 --ie_filter 0.5 --num_reps 8`
#predict
`python /private/groups/brookslab/althornt/apps/eVIP/bin/eVIP_predict.py -i /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/zscore_eVIP/eVIP_compare_out.txt -o /private/groups/brookslab/althornt/RNF43-eVIP-analysis/L1000-wtcs-vs-spearman/zscore_eVIP/eVIP_predict_out --mut_wt_rep_thresh 0.05 --mut_wt_rep_rank_diff 0 --disting_thresh 0.05 --use_c_pval`

"""

#open eVIP predict results from wtcs
wtcs_predict = pd.read_csv("wtcs_eVIP_predict_out.txt",sep="\t",index_col='mut',engine='python')
wtcs_predict_df = wtcs_predict['prediction']
wtcs_predict_pval = wtcs_predict['wt_mut_rep_vs_wt_mut_conn_c_pval']
wtcs_predict_pval = wtcs_predict_pval.rename(columns={"wt_mut_rep_vs_wt_mut_conn_c_pval":"wtcs"})

#import eVIP paper results
paper = pd.read_csv("eVIP_paper_L1000_calls.txt",sep="\t",index_col="Allele")
paper_df = paper['eVIP_prediction']

#compare paper results to wtcs results
paper_vs_wtcs = pd.merge(paper_df.to_frame(), wtcs_predict_df.to_frame(), left_index=True, right_index=True)
print(paper_vs_wtcs.head())

#open predict results from z-score, spearman
zscore_predict = pd.read_csv("spearman_eVIP_predict_out.txt",sep="\t",index_col='mut')
zscore_predict_df = zscore_predict['prediction']
zscore_predict_pval = zscore_predict['wt_mut_rep_vs_wt_mut_conn_c_pval']
zscore_predict_pval = zscore_predict_pval.rename(columns={"wt_mut_rep_vs_wt_mut_conn_c_pval":"zscore"})

#compare paper results  to zscore,  spearman
paper_vs_zscore = pd.merge(paper_df.to_frame(), zscore_predict_df.to_frame(), left_index=True, right_index=True)
print(paper_vs_zscore.head())

#prediction zscore vs wtcs
wtcs_vs_zscore_prediction = pd.merge(wtcs_predict_df.to_frame(), zscore_predict_df.to_frame(), left_index=True, right_index=True)
wtcs_vs_zscore_prediction[wtcs_vs_zscore_prediction['prediction_x']==wtcs_vs_zscore_prediction['prediction_y']].shape

#set up df for plot
impactful = ['COF','GOF','LOF']
notimpactful = ['NI','Neutral']

for i in wtcs_vs_zscore_prediction.index:
    if wtcs_vs_zscore_prediction.loc[i,'prediction_x'] in impactful and wtcs_vs_zscore_prediction.loc[i,'prediction_y'] in impactful:
        wtcs_vs_zscore_prediction.at[i, 'consistent prediction'] = "yes"

    elif wtcs_vs_zscore_prediction.loc[i,'prediction_x'] in notimpactful and wtcs_vs_zscore_prediction.loc[i,'prediction_y'] in notimpactful:
        wtcs_vs_zscore_prediction.at[i, 'consistent prediction'] = "yes"

    else:
        wtcs_vs_zscore_prediction.at[i, 'consistent prediction'] = " no"

print(wtcs_vs_zscore_prediction.head())

wtcs_vs_zscore_pval = pd.merge(wtcs_predict_pval.to_frame(), zscore_predict_pval.to_frame(), left_index=True, right_index=True)
wtcs_vs_zscore_color = pd.merge(-np.log10(wtcs_vs_zscore_pval), wtcs_vs_zscore_prediction, left_index=True, right_index=True)

#Use lmplot to plot scatter points
graph = sns.lmplot(x='0_x', y='0_y', hue='consistent prediction', data=wtcs_vs_zscore_color, fit_reg=False, scatter_kws={'alpha':0.55})
#Use regplot to plot the regression line for the whole points
sns.regplot(x='0_x', y='0_y', data=wtcs_vs_zscore_color, scatter=False, ax=graph.axes[0, 0])
plt.title('Adjusted p-values from L1000 wtcs vs spearman')
plt.xlabel('-log10(pval) from wtcs')
plt.ylabel('-log10(pval) from spearman')
#add vertical line for significance
plt.axvline(x=-np.log10(.05), color = 'grey', linestyle = '--')
#calculating pearson and rounding
r2 = str(stats.pearsonr(wtcs_vs_zscore_color['0_x'], wtcs_vs_zscore_color['0_y'])[0])[:5]
#box style
props = dict(boxstyle='square', facecolor='grey', alpha=0.1)

plt.text(1.5, 4, "pearson="+r2, fontsize=12,bbox=props)

plt.savefig("wtcs-vs-spearman.svg", dpi = 400, bbox_inches = "tight")
plt.clf()
plt.close()
