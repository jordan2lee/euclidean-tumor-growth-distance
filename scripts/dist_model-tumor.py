# Calculate Distance Between Tumor and Derived Model
# using Custom OHSU Euclidean Distance Method

import pandas as pd
import numpy as np
import statistics
import scipy.stats as stats
import json
from collections import Counter

#######
cancer = 'PAAD'
run = 'HCMITumorModel'
matches_file = 'src/distance_metric/HCMI_AWG_Model-Tumor-Normal_Linkage_v2.0_2.20.2024.txt'
ensure_matches_file = 'data/distance_metric/src/missing_ohsu_euclidean_distance_all.02.12.2024.csv'
#######

def df_read(CANCER, SAMPLE_TYPE):
    '''
    read in appropriate file and filter for TMP 
    important genes. provide HCMI abbrev for cancer cohort
    and sample type (Model or Tumor).
    outputs filtered df
    '''
    # if needed, merge the two matrices 
    if CANCER == 'LUNG' or CANCER == 'ESO':
        if CANCER == 'LUNG':
            file = 'data/midway.freeze.v2/{}_GEXP/{}_GEXP_prep2_{}.tsv'.format(CANCER, 'LUAD', SAMPLE_TYPE)
            df = pd.read_csv(file, sep='\t', index_col=0).sort_index(axis=0)
            file = 'data/midway.freeze.v2/{}_GEXP/{}_GEXP_prep2_{}.tsv'.format(CANCER, 'LUSC', SAMPLE_TYPE)
            df2 = pd.read_csv(file, sep='\t', index_col=0).sort_index(axis=0)
            assert list(df.index)==list(df2.index) # req for df concat using
        elif CANCER == 'ESO':
            file = 'data/midway.freeze.v2/{}_GEXP/{}_GEXP_prep2_{}.tsv'.format(CANCER, 'ESCC', SAMPLE_TYPE)
            df = pd.read_csv(file, sep='\t', index_col=0).sort_index(axis=0)
            file = 'data/midway.freeze.v2/{}_GEXP/{}_GEXP_prep2_{}.tsv'.format(CANCER, 'GEA', SAMPLE_TYPE)
            df2 = pd.read_csv(file, sep='\t', index_col=0).sort_index(axis=0)
            assert list(df.index)==list(df2.index) # req for df concat using        

        # Add any new fts from second CANCER cohort that isn't in first cohort
        b= list(df2.columns)
        a= list(df.columns)
        new = [ft for ft in b if ft not in a]
        if len(new) !=0:
            s1 = df2[new]
            df =pd.concat([df, s1], join="outer", axis = 1)
    else:
        file = 'data/midway.freeze.v2/{}_GEXP/{}_GEXP_prep2_{}.tsv'.format(CANCER, CANCER, SAMPLE_TYPE)
        df = pd.read_csv(file, sep='\t', index_col=0)
        if CANCER == 'LGGGBM':
            if 'HCM-BROD-1124-C16-06A' in list(df.index):
                print('applying manual removal for HCM-BROD-1124-C16-06A in LGGGBM')
                df = df.drop(['HCM-BROD-1124-C16-06A' ], axis=0, )
            if 'HCM-BROD-1124-C16-85B' in list(df.index):
                print('applying manual removal for HCM-BROD-1124-C16-85B in LGGGBM')
                df = df.drop(['HCM-BROD-1124-C16-85B' ], axis=0, )
    return df


def intersection(list1, list2):
    overlap = [a for a in list1 if a in list2]
    return overlap



def get_euclidean(metric, v1, v2):
    '''
    example: ('single', 1, 6) indicates single euclidean distance between value 1 and value 2    
    example: ('mean', list(1,3,4,5), list(6,7,8,9)) indicates euclidean distance between each pair of values (order matters)
    and then the mean euclidean distance of those computed
    '''
    import numpy as np
    import statistics
    # calcuate the eucl distance for pair of values
    if metric == 'single':
        euc = np.linalg.norm(v1 - v2)
        return round(euc, 2)
    
    # calcuate the mean eucl distance for the pair of lists
    if metric == 'mean':
        assert len(v1) == len(v2), 'two lists must be same length'
        euc_list = []
        for i in range(0, len(v1)):
            euc = np.linalg.norm(v1[i] - v2[i])
            euc_list.append(euc)
        mean_euc = statistics.mean(euc_list) 
        return round(mean_euc, 2)

# -----
# Prep
# -----

# Get important feature lists (important genes = TMP ft selected sets)
with open('data/distance_metric/src/cancer2fts.json', 'r') as fh:
    for line in fh:
        cancer2fts = json.loads(line)

# read in hcmi data
model_df = df_read(cancer, 'Model')
tumor_df = df_read(cancer, 'Tumor')

# ft reduction (using TMP genes) and filter for shared genes (ex. tumor, model, TMP impt genes) 
genes_keep = intersection( list(tumor_df.columns), list(model_df.columns) )
combined_top_fts = cancer2fts[cancer] # TMP important fts
genes_keep = intersection( genes_keep, combined_top_fts)
model_df = model_df[genes_keep]
tumor_df = tumor_df[genes_keep]
# now save these files
model_df.to_csv('data/distance_metric/{}_GEXP/{}_topgenes_Model.tsv'.format(cancer, cancer), sep='\t', index=True)
tumor_df.to_csv('data/distance_metric/{}_GEXP/{}_topgenes_Tumor.tsv'.format(cancer,cancer), sep='\t', index=True)

# -----
# Calculate distance
# -----

### Distance between a given tumor and the rest of it's tumor cohort ###
### ex. What is the distance between PAAD_tumor_1 and the rest of the PAAD_tumor cohort? ###
# Calculate the overall mean euc_dist sample 1 is to all other samples (find mean dist in table)
results_sample = []
results_distance = []
results_type = []
results_sample_info = []
for sample_1 in tumor_df.index:
    samples = list(tumor_df.index)
    samples.remove(sample_1) # ensure not comparing the same sample

    # Compute mean sample distance of ONE sample to each other sample (pairwise)
    list_sample_1 =[]
    list_sample_2 =[]
    euclidean_list =[]
    for sample_2 in samples:
        # get expression values across all genes
        value_1 = list(tumor_df.loc[sample_1])
        value_2 = list(tumor_df.loc[sample_2])
        genes = list(tumor_df.columns)

        # calculate mean euc distance between 2 samples (mean across all genes)
        euc_mean = get_euclidean('mean', value_1, value_2)

        # prep for results matrix
        list_sample_1.append(sample_1)
        list_sample_2.append(sample_2)
        euclidean_list.append(euc_mean)

    # View results
    # pd.DataFrame(list(zip(list_sample_1, list_sample_2, euclidean_list)), columns = ['sample1', 'sample2', 'mean_euclidean_dist'])

    # Calculate the overall mean euc_dist sample 1 is to all other samples (find mean dist in table)
    mean_euc_1vAll_samples = round(statistics.mean(euclidean_list), 2)

    # Save
    results_sample.append(sample_1)
    results_distance.append(mean_euc_1vAll_samples)
    results_type.append('mean_euc_1vsALL_samples')
    results_sample_info.append('tumorVStumors_sameTissue')
results_df = pd.DataFrame(list(zip(results_sample, results_distance, results_type, results_sample_info)), 
             columns = ['sample_1', 'distance', 'distance_type', 'sample_info'])
# calculate z-score (within tumors )
darray = np.array(results_df['distance'])
z = stats.zscore(darray)
results_df['z-score_to_tumors']= z

### Distance between a given tumor and it's derived model ###
### ex. What is the distance between PAAD_tumor_1 and the PAAD_model_1 ? ###

# Create look up matches of model and cancer model
ref = pd.read_csv(matches_file, sep='\t')

# note some tumors match to multiple models
tumor2model = ref[['Matched Tumor Aliquot', 'Cancer Model Aliquot']]
tumor2model = tumor2model.dropna() #only keep tum:model pairs that have no nan

# tweak based on missing file 
# add any pairs that are missing
mdf = pd.read_csv(ensure_matches_file, sep=',')

# Create abbrev to full cancer cohort dict
with open('src/full_cancers.json', 'r') as fh:
    for line in fh:
        longform = json.loads(line)
        
mdf= mdf[mdf['cancer_type'].isin(longform[cancer])]
mdf = mdf[['matched_tumor_aliquot', 'aliquot_id3', ]].reset_index(drop=True)
mdf.columns = ['Matched Tumor Aliquot','Cancer Model Aliquot']

# Build on "matched tumor aliquot" column from Missing values file
new = {'Matched Tumor Aliquot':[], 'Cancer Model Aliquot':[]}
for i in range(0, mdf.shape[0]):
    mdf_tumor = mdf['Matched Tumor Aliquot'][i]
    mdf_model = mdf['Cancer Model Aliquot'][i]
    # if tumor not in tracker, then add it
    if mdf_tumor not in list(tumor2model['Matched Tumor Aliquot']):
        new['Matched Tumor Aliquot'].append(mdf_tumor)
        new['Cancer Model Aliquot'].append(mdf_model)
    # if tumor has a different model match, then add it
    elif list(tumor2model[tumor2model['Matched Tumor Aliquot']==mdf_tumor]['Cancer Model Aliquot'])[0] != mdf_model:
        new['Matched Tumor Aliquot'].append(mdf_tumor)
        new['Cancer Model Aliquot'].append(mdf_model)
new_df = pd.DataFrame.from_dict(new)
tumor2model = pd.concat([tumor2model, new_df])
tumor2model = tumor2model.reset_index(drop=True)

# Build on "matched tumor aliquot" column from Missing values file
new = {'Matched Tumor Aliquot':[], 'Cancer Model Aliquot':[]}
for i in range(0, mdf.shape[0]):
    mdf_tumor = mdf['Matched Tumor Aliquot'][i]
    mdf_model = mdf['Cancer Model Aliquot'][i]
    # if tumor not in tracker, then add it
    if mdf_model not in list(tumor2model['Cancer Model Aliquot']):
        new['Matched Tumor Aliquot'].append(mdf_tumor)
        new['Cancer Model Aliquot'].append(mdf_model)
        print('adding new model {}'.format(mdf_model))
    # if tumor has a different model match, then add it
    elif list(tumor2model[tumor2model['Cancer Model Aliquot']==mdf_model]['Matched Tumor Aliquot'])[0] != mdf_tumor:
        new['Matched Tumor Aliquot'].append(mdf_tumor)
        new['Cancer Model Aliquot'].append(mdf_model)
        print('model {} will now have tumor match {} instead of {}'.format(mdf_model, mdf_tumor,list(tumor2model[tumor2model['Cancer Model Aliquot']==mdf_model]['Matched Tumor Aliquot'])[0] ))
new_df = pd.DataFrame.from_dict(new)
tumor2model = pd.concat([tumor2model, new_df])
tumor2model = tumor2model.reset_index(drop=True)

# Create tumor-model mappings that need to be forced  
with open('src/extra_tumor_model_pairs.json', 'r') as fh:
    for line in fh:
        manual_matches = json.loads(line)

# Add this if not already there or if already there then update the pairings
for ipair in manual_matches[cancer]:
    sample_tumor = ipair['tumor']
    sample_model = ipair['model']
    if sample_tumor not in list(tumor2model['Matched Tumor Aliquot']):
        tumor2model.loc[len(tumor2model.index)] = [sample_tumor, sample_model]
        
    if sample_model not in list(tumor2model['Cancer Model Aliquot']):
        tumor2model.loc[len(tumor2model.index)] = [sample_tumor, sample_model]

list_sample_1 =[]
list_sample_2 =[]
euclidean_list =[]

for sample in results_df['sample_1']: # results_df is of HCMI tumors
    model_from_tumor2model = list(tumor2model[tumor2model['Matched Tumor Aliquot']==sample]['Cancer Model Aliquot'])
    # loop thru each pair (if tumor matches to multiple models)
    for i_model in model_from_tumor2model:
        if i_model in list(model_df.index):
            # get expression values across all genes
            value_1 = list(tumor_df.loc[sample, ])
            value_2 = list(model_df.loc[i_model ,])
            genes = list(tumor_df.columns)
            
            # calculate mean euc distance between 2 samples (mean across all genes)
            euc_mean = get_euclidean('mean', value_1, value_2)

            # prep for results matrix
            list_sample_1.append(sample)
            list_sample_2.append(i_model)
            euclidean_list.append(euc_mean)

            ## View results
            #pd.DataFrame(list(zip(list_sample_1, list_sample_2, euclidean_list)), columns = ['sample1', 'sample2', 'mean_euclidean_dist'])


results_df_2 = pd.DataFrame(list(zip(list_sample_1, list_sample_2, euclidean_list, ['mean_euc']*len(list_sample_1), ['tumorVSmodel_sameTissue']*len(list_sample_1) )), 
             columns = ['sample_1', 'sample_2', 'distance', 'distance_type', 'sample_info'])

# calculate z-score (within tumors )
darray = np.array(results_df_2['distance'])
z = stats.zscore(darray)
results_df_2['z-score_to_model']= z

# merge 2 df to one easy to ready matrix
# drops tumor row if not found in both df
a = results_df[['sample_1','distance','sample_info', 'z-score_to_tumors']]
a.columns = ['tumor', 'dist_1', 'dist_1_info','z-score_to_tumors']
b = results_df_2[['sample_1', 'sample_2', 'distance', 'sample_info', 'z-score_to_model']]
b.columns = ['tumor', 'model', 'dist_2', 'dist_2_info','z-score_to_model']
final_results = pd.merge(a,b, on='tumor', how='inner')

# now compare dist of model to other_tumors_in_cohort
# <1 model is more similar to tumor (than tumor to tumor)
# 1 tumor and model are the same
# >1 model is less similar to tumor (than tumor to tumor)
final_results['ratio'] = final_results['dist_2']/final_results['dist_1']
final_results = final_results[['tumor', 'model', 'ratio', 'dist_1', 'dist_1_info', 'dist_2', 'dist_2_info','z-score_to_tumors','z-score_to_model']]

# determine which are similar based on z-score
z_tumors_list = list(final_results['z-score_to_tumors'])
z_models_list = list(final_results['z-score_to_model'])
similar = []
for i in range(0, len(z_tumors_list)):
    z_tum = z_tumors_list[i]
    z_mod = z_models_list[i]
    if abs(z_tum) <= 3:
        if abs(z_mod) <= 3:
            similar.append('similarTumor:similarModel')
        else:
            similar.append('similarTumor:dissimilarModel')
    else:
        if abs(z_mod) <= 3:
            similar.append('dissimilarTumor:similarModel')
        else:
            similar.append('dissimilarTumor:dissimilarModel')
            
final_results['similar']= similar
final_results.to_csv('data/distance_metric/main_results/euc_ratio_HCMITumor.Model_{}.tsv'.format(cancer), sep='\t', index=False)

# Create outlier file
problem_labels = ['similarTumor:dissimilarModel','dissimilarTumor:similarModel', 'dissimilarTumor:dissimilarModel']
issues =final_results[final_results['similar'].isin(problem_labels)].reset_index(drop=True)
issues.to_csv('data/distance_metric/main_results/outlier_samples_HCMITumor.Model_{}.tsv'.format(cancer), sep='\t',index=False)