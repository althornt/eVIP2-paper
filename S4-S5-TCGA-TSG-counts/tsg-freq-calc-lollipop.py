import pandas as pd
from collections import Counter
import os

"""
files from cbioportal RNF43 mutation downloads

input TUSON_TSG_OG_selected_predictions_141007_TSG_q0.05.txt into cbioportal to get mutations for each sample for each mutations
had to do it in batches
only TCGA pancancer studies
"""

def  main():
    #import different batches of TSG cbioportal mutation files
    batch1 = pd.read_csv("input_data/cbioportal_TSG_mut_batch1.txt", sep="\t",index_col ="SAMPLE_ID")
    batch2 = pd.read_csv("input_data/cbioportal_TSG_mut_batch2.txt", sep="\t",index_col ="SAMPLE_ID")
    batch3 = pd.read_csv("input_data/cbioportal_TSG_mut_batch3.txt", sep="\t",index_col ="SAMPLE_ID")

    #combine batches
    TSG_mut_df = pd.concat([batch1,batch2,batch3], axis=1)
    TSG_mut_df = TSG_mut_df.loc[:,~TSG_mut_df.columns.duplicated()]

    ###############################################################################
    #per cohort mutation percent calculation

    output_list = calc(TSG_mut_df)

    #list of lists to df
    FS_df_combine = pd.DataFrame(output_list)
    FS_df_combine.columns = ['gene','tcga_cohort','mutation','percent_in_cohort','count_in_cohort']
    FS_df_combine = FS_df_combine.sort_values(by="gene",ascending=True)
    FS_df_combine.to_csv("TSG_tcga_cohort_freq_by_gene.csv")
    FS_df_combine = FS_df_combine.sort_values(by="count_in_cohort",ascending=False)
    FS_df_combine.to_csv("TSG_tcga_cohort_freq_by_count.csv")

    #################
    # genes over 5%
    #################
    print("Over 5% frequency")
    results_over_5percent_1count = FS_df_combine[(FS_df_combine["percent_in_cohort"]>5) & (FS_df_combine["count_in_cohort"]>1)]
    print(results_over_5percent_1count)

    gene2cohort = []
    #group results by gene for plotting only the cohorts
    groupbygene = results_over_5percent_1count.groupby('gene')
    for key, item in groupbygene:
        print(key)
        muts = set(groupbygene.get_group(key)["mutation"].tolist())

        cohorts = groupbygene.get_group(key)["tcga_cohort"].tolist()
        gene2cohort.append([key,muts,cohorts])

    lolli(TSG_mut_df,gene2cohort)


    #################
    # 2 genes over 1%
    #################
    gene2cohort_2 = []

    results_over_1percent_1count = FS_df_combine[(FS_df_combine["percent_in_cohort"]>1) & (FS_df_combine["count_in_cohort"]>1)]
    print(results_over_1percent_1count)
    result_group = results_over_1percent_1count.groupby("gene")

    for key, item in result_group:
        if len(set(result_group.get_group(key)["mutation"].tolist())) > 1 :
            print(key)

            muts = set(result_group.get_group(key)["mutation"].tolist())
            cohorts = set(result_group.get_group(key)["tcga_cohort"].tolist())
            gene2cohort_2.append([key,muts,cohorts])

    print(gene2cohort_2)
    lolli(TSG_mut_df,gene2cohort_2)


def lolli(TSG_mut_df,gene2cohort):
    for i in gene2cohort:
        gene = i[0]
        muts = [i.split("*")[0].split("fs")[0][:-1]+"fs" for i in i[1]]
        print(muts)

        cohort = i[2]

        df = TSG_mut_df[TSG_mut_df["STUDY_ID"].isin(cohort)]
        df = df[[gene]]


        l = df[gene].tolist()
        l = [i for i in l  if str(i) != 'nan']

        l = [i.split() for i in l]
        flat_list = [item for sublist in l for item in sublist if item not in ["WT","NS"] ]
        li = [i.split("*")[0] for i in flat_list ]

        #remove  extra character before fs so in order to combine fs the same position
        new_list = []

        #only plot frame shift , rename so frameshifts at the same position are combined
        for  i in li:
            if i.endswith("fs"):
                new_list.append(i.split("fs")[0][:-1]+"fs")

            # else:
            #     new_list.append(i)

        count  = Counter(new_list)
        print(count)

        cmd = "./lollipops-v1.5.2-mac64/lollipops  -domain-labels=off    -w=1000 -o ./lollipop_out/"+gene+"_out.svg  "+gene+" "

        for key, val in count.items():
            if val == 1 :
                size=0.1
            else:
                size=val**2

            if key  in muts:
                color= "#0000ff"

            elif key.split("*")[0].endswith("fs"):
                color= "#228B22"

            else:
                color = "#000000"

            cmd += str(key)+color+"@"+str(size)+" "

        print(cmd)
        os.system(cmd)



def calc(TSG_mut_df):
    output_list = []

    for gene in TSG_mut_df.columns:
        #group by cancer type
        grouped = TSG_mut_df.groupby('STUDY_ID')

        #for each cancer type
        for name,group in grouped[gene]:
            samples_in_cohort = float(len(group))
            # print(name,samples_in_cohort)

            #count mutation occurances
            df = group.value_counts(dropna=True).to_frame()

            #if there is more than one mutation add count to the single mutation
            for i  in df.index:
                if len(i.split()) > 1:
                    for mut in i.split():
                        try:
                            df.loc[mut] +=1
                        except:
                            df.loc[mut] =1
                    #then remove the row with multiple mutations
                    df = df.drop([i])

            for x  in df.index:
                if "fs" in x:
                    percent =  (float(df.loc[x])/samples_in_cohort)*100
                    if percent >1:
                        # print(gene,name,x,percent,int(df.loc[x]))
                        output_list.append([gene,name,x,percent,int(df.loc[x])])

                        # print(FS_df_combine)
    return output_list

if __name__ == '__main__':
    main()
