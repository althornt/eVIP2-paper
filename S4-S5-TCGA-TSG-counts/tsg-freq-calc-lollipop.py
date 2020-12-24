import pandas as pd
from collections import Counter
import os

"""
> cbioportal.com
> Select all TCGA and DFCI cohorts
> Input list of TSG into cbioportal "TUSON_TSG_OG_selected_predictions_141007_TSG_q0.txt"
> There is a limit, so have to do it in batches
> Download mutation file for each batch
"""

def  main():
    #import different batches of TSG cbioportal mutation files
    batch1 = pd.read_csv("input_data/TCGA_DFCI/cbioportal_TSG_mut_TCGA_DFCI_batch1.txt", sep="\t",index_col ="SAMPLE_ID")
    batch2 = pd.read_csv("input_data/TCGA_DFCI/cbioportal_TSG_mut_TCGA_DFCI_batch2.txt", sep="\t",index_col ="SAMPLE_ID")
    batch3 = pd.read_csv("input_data/TCGA_DFCI/cbioportal_TSG_mut_TCGA_DFCI_batch3.txt", sep="\t",index_col ="SAMPLE_ID")

    """
    batch1 = pd.read_csv("input_data/cbioportal_TSG_mut_batch1.txt", sep="\t",index_col ="SAMPLE_ID")
    batch2 = pd.read_csv("input_data/cbioportal_TSG_mut_batch2.txt", sep="\t",index_col ="SAMPLE_ID")
    batch3 = pd.read_csv("input_data/cbioportal_TSG_mut_batch3.txt", sep="\t",index_col ="SAMPLE_ID")
    """

    #combine batches
    TSG_mut_df = pd.concat([batch1,batch2,batch3], axis=1)
    TSG_mut_df = TSG_mut_df.loc[:,~TSG_mut_df.columns.duplicated()]

    # print(set(TSG_mut_df["STUDY_ID"].tolist()))

    #per cohort mutation percent calculation
    freq_calc_output_lists = cohortCalc(TSG_mut_df)
    # print(freq_count_dfs)

    #Combining results - list of lists to df
    FS_df_combine = pd.DataFrame(freq_calc_output_lists)
    FS_df_combine.columns = ['gene','tcga_cohort','mutation','percent_in_cohort','count_in_cohort','df']
    FS_df_combine = FS_df_combine.sort_values(by="gene",ascending=True)
    FS_df_combine.drop(["df"],axis=1).to_csv("TSG_tcga_cohort_freq_by_gene.csv" ,index=False)
    FS_df_combine = FS_df_combine.sort_values(by="count_in_cohort",ascending=False)
    FS_df_combine.drop(["df"],axis=1).to_csv("TSG_tcga_cohort_freq_by_count.csv", index=False)


    #################
    # genes over 5%
    #################
    print("Over 5% frequency")
    results_over_5percent_1count = FS_df_combine[(FS_df_combine["percent_in_cohort"]>5) & (FS_df_combine["count_in_cohort"]>1)]
    print(results_over_5percent_1count.drop(["df"],axis=1))

    gene2cohort = []
    #group results by gene for plotting only the cohorts
    groupbygene = results_over_5percent_1count.groupby('gene')
    for key, item in groupbygene:
        muts = list(set(groupbygene.get_group(key)["mutation"].tolist()))
        cohorts = groupbygene.get_group(key)["tcga_cohort"].tolist()
        dfs = groupbygene.get_group(key)["df"].tolist()

        gene2cohort.append([key,muts,cohorts,dfs])

    # print(groupbygene)
    # for i in gene2cohort:
    #     print(i[0:3])
    lolli(gene2cohort,"lollipops-FS-freq-over5")


    ###########################
    # genes with 2 FS over 1%
    ###########################
    gene2cohort_2 = []

    results_over_1percent_1count = FS_df_combine[(FS_df_combine["percent_in_cohort"]>1) & (FS_df_combine["count_in_cohort"]>1)]
    # print(results_over_1percent_1count)
    result_group = results_over_1percent_1count.groupby("gene")

    for key, item in result_group:
        if len(set(result_group.get_group(key)["mutation"].tolist())) > 1 :
            muts = list(set(result_group.get_group(key)["mutation"].tolist()))
            cohorts = list(set(result_group.get_group(key)["tcga_cohort"].tolist()))
            dfs = result_group.get_group(key)["df"].tolist()

            #get df that for unique cohorts for plots
            dfs_clear = result_group.get_group(key).drop_duplicates('tcga_cohort')["df"].tolist()

            gene2cohort_2.append([key,muts,cohorts, dfs_clear])


    # for i in gene2cohort_2:
    #     print(i[0:2])

    lolli(gene2cohort_2, "lollipops-2-FS-freq-over1")



def lolli(gene2cohort,outdir):
    for i in gene2cohort:
        gene,muts,cohort,dfs = i[0],i[1],i[2],i[3]

        #combine dfs from different cohorts
        comb_df =  pd.concat(dfs, axis=1)
        comb_df["sum"] =comb_df.sum(axis=1)
        count_dict = comb_df["sum"].to_dict()

        cmd = "./lollipops-v1.5.2-mac64/lollipops  -domain-labels=off    -w=1000 -o ./"+outdir+"/"+gene+"_out.svg  "+gene+" "

        for key, val in count_dict.items():
            if val == 1 :
            #     size=0.1
            # else:
            #     size=val**2

            size=val

            if key  in muts:
                color= "#0000ff"

            elif key.split("*")[0].endswith("fs"):
                color= "#228B22"

            else:
                color = "#000000"

            cmd += str(key)+color+"@"+str(size)+" "

        print(cmd)
        os.system(cmd)


def cohortCalc(TSG_mut_df):
    """
    Use df of mutations per sample to calculate frequency of each FS mutation
    """

    freq_calc_output = []
    for gene in TSG_mut_df.columns:
        #group by cancer type
        grouped = TSG_mut_df.groupby('STUDY_ID')

        #for each cancer cohort
        for name,group in grouped[gene]:
            samples_in_cohort = float(len(group))

            #count mutation occurances
            df = group.value_counts(dropna=True).to_frame()
            df_fs_only =  pd.DataFrame(columns = df.columns)

            #for each mut , if there is more than one mutation add count to the single mutation
            for i in df.index:
                #only 1 mutation , rename
                if len(i.split()) ==1:
                    if "fs*" in i:
                        fs_mut =i.split("fs")[0][:-1]+"fs"
                        try:
                            df_fs_only.loc[fs_mut] += df.loc[i]
                        except:
                            df_fs_only.loc[fs_mut] = df.loc[i]

                #if more than 1 mutation
                elif len(i.split()) > 1:
                    for mut in i.split():

                        #if frameshift, rename (R117Afs*41 and R117Pfs*41 --> R117fs)
                        if "fs*" in mut:
                            fs_mut =mut.split("fs")[0][:-1]+"fs"

                            try:
                                df_fs_only.loc[fs_mut] += df.loc[i]
                            except:
                                df_fs_only.loc[fs_mut] = df.loc[i]

            #calculate frequency of each FS mutation in each cohort
            for x  in df_fs_only.index:
                if "fs" in x:
                    percent =  (float(df_fs_only.loc[x])/samples_in_cohort)*100
                    freq_calc_output.append([gene,name,x,percent,int(df_fs_only.loc[x]),df_fs_only])

    return freq_calc_output


if __name__ == '__main__':
    main()
