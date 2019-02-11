from lifelines.datasets import load_rossi
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
import os
from os import path
from glob import glob
import numpy as np
import pandas as pd
from tqdm import tqdm, trange
from utils import *
import sys, io
import warnings

model = 'mutex_t07_nsep_cov'

out_path = '../hint/out/'
data_path = '../data/'
clinic_path = data_path+'clinical_data/'
clinic_file = '{}_clinical_data.txt'
exp_gene_file = '../data/TCGA/{}_Normalized.csv'
# List of cancers studied
cancer_types = ['BLCA','BRCA','CRC', 'GBM', 'HNSC', 'KIRC', 'LAML', 'LUAD', 'LUSC', 'OV', 'UCEC']

with open('../data/exp_final_gene_list.txt', 'r') as f:
    final_genes= [l.rstrip() for l in f.readlines()]
with open('../data/hint_inters_pan.txt', 'r') as f:
    hint_u_pan= [l.rstrip() for l in f.readlines()]

comps = []
significant_cancers =[]
with open('../hint/out/cancer_subtype_test_results/{}/cc_n100_cancer_subtype_tests.txt'.format(model), 'r') as f:
    lines = f.readlines()[:]
    for i in range(0,len(lines), 12):
        comp = lines[i].strip().split('[')[-1].split(',')
        comp = [c.replace(']', '').replace('\'', '').strip() for c in comp]
        comps.append(comp)
        type_interest = [j-1 for j in range(1,12) if len(lines[i+j].split()) ==8]
        significant_cancers.append(type_interest)

# This is a fix for 'MLL2', 'MLL3', 'PHF17' which are not in the mapping file
comps[7][1:] = ['KMT2B', 'KMT2C']
comps[-1][8] = 'JADE1'

cancer2modules = {}
for i,m in enumerate(comps):
    for t in significant_cancers[i]:
        if cancer_types[t] in cancer2modules.keys():
            cancer2modules[cancer_types[t]].append(m)
        else:
            cancer2modules[cancer_types[t]] = [m]

module2id = {str(m):i for i,m in enumerate(comps)}
id2module = {i:str(m) for i,m in enumerate(comps)}




genes = []
for c in comps:
    genes.extend(c)
# df_m: pandas Dataframe contining the mapping between the gene id and gene symbol
df_m = pd.read_csv('../data/ensembl_biomart_hg19_proteincoding_gene_refseq_mapping.txt',delimiter = '\t' )

df_m = df_m[['Gene stable ID','Gene name' ]]

# df_m1: pandas Dataframe where the Gene name column is the index
df_m1 = df_m[df_m['Gene name'].isin(genes)].drop_duplicates('Gene stable ID')
df_m1 = df_m1.set_index('Gene name')


# df_m2: pandas Dataframe where the Gene ID column is the index
df_m2 = df_m1.reset_index().set_index('Gene stable ID')

# gene_ids: contins the gene ids that map the genes of interest (genes in our modules) to their symbols
gene_ids = df_m2.index.values


# Loop over all the cancers, and run the analysis on the modules where cancer_ is significant

for cancer_ in cancer_types:
    # the clinical file for CRC cancer starts with 'COADREAD' instead of CRC
    if cancer_ == 'CRC':
        clinical_data = clinic_file.format('COADREAD')
    else:
        clinical_data = clinic_file.format(cancer_)

    # CLINICAL DATA

    # Load the clinical data for the cancer_
    df_clinical = pd.read_csv(clinic_path+clinical_data,usecols = ['patient_id','event','time_to_event'], delimiter = '\t', low_memory=False)

    df_clinical = df_clinical.set_index('patient_id')
    df_clinical.index.names = ['patients']

    # Safe check: Delete any rows with nan entries
    df_clinical =df_clinical[~df_clinical.isnull().any(axis=1)]

    # GENE EXPRESSION DATA
    # laod the gene expression for the cancer type. (The data was doawnloaded using R code)
    # The data contains the normalized raw counts
    if cancer_ == 'CRC':
        df = pd.read_csv(exp_gene_file.format('COAD'), low_memory=False)
    else:
        df = pd.read_csv(exp_gene_file.format(cancer_), low_memory=False)

    df = df.set_index('Unnamed: 0')

    # select only the genes of interest
    df_s = df.loc[gene_ids]

    # concatenate the mapping data with gene expression
    # the join is done on the index (gene id)
    df_s = df_s.join(df_m2)

    # set the gene name as the index
    df_s = df_s.set_index('Gene name')

    # group the rows by gene name, and take the mean (this is due to multiple mapping from gene id to gene name)
    # but for each dataset there is only one valid mapping in the gene expression,
    # so the mean serves as away to select the valid entry
    df_s_group = df_s.groupby(level=0)
    df_s_group = df_s_group.mean()

    # Transpose the matrix (each row is patient, and each column is a gene)
    df_s_group = df_s_group.T

    df_s_group.index.names = ['patients']
    df_s_group.columns.name = 'genes'

    # Select the patients with tumor
    # For LAML (blood cancer), the tumor patients have a tumor ID of 03X (X= A,B,...)
    if cancer_ == 'LAML':
        keep_patient = [x.split('-')[3][:2] == '03' for x in df_s_group.index]
    else:
        keep_patient = [x.split('-')[3][:2] == '01' for x in df_s_group.index]

    # select the patients
    df_s_group = df_s_group.loc[keep_patient]

    # Rename the patients using only the first three identifiers
    df_s_group.index = ['-'.join(x.split('-')[:3]) for x in df_s_group.index]
    df_s_group.index.names = ['patients']

    # Concatenate the gene expression data with the clinical data
    # join is done on the patients id
    df2 = df_clinical.join(df_s_group)
    df2 = df2.loc[~df2[df2.columns[2:]].isnull().any(axis=1)]
    df2 = df2.groupby(df2.index).first()

    # for each module that is significant in the given cancer_, do the analysis
    for module in cancer2modules[cancer_]:

        risk_summary = []

        # m_id is the id given for each module, based on the 'cc_n100_cancer_subtype_tests.txt' file
        m_id = module2id[str(module)]

        summaries = []
        results = []
        gene2coef = {}
        # Calculate the Coxph coeff for each gene in the module
        for g in module:

            tmp_df  = df2[[g]+['event','time_to_event']].astype('float')
            cph = CoxPHFitter()
            cph.fit(tmp_df, duration_col='time_to_event', event_col='event', show_progress=False)

            stdout = sys.stdout
            sys.stdout = io.StringIO()
            cph.print_summary(decimals=5)
            output = sys.stdout.getvalue()
            sys.stdout = stdout
            gene2coef[g] = cph.summary[['coef']].values
            summaries.append(output+'\n\n')

        if not os.path.exists('../results/CoxPh/{}/'.format(m_id)):
                os.mkdir('../results/CoxPh/{}/'.format(m_id))
        # Write the coefficients to file
        # with open('../results/CoxPh/{}/CoxPh_{}_summary.txt'.format(m_id, cancer_), 'w') as f:
        #     f.write(''.join(summaries))

        # Risk Analysis

        # Split the data into 50% training and 50% testing
        train_idx = df2.sample(frac= 0.5, random_state=1234).index

        df_train = df2.loc[train_idx]
        df_test = df2.loc[~df2.index.isin(train_idx)]

        coiffiecients = [gene2coef[f] for f in module]
        # Calculate the S-score for training and testing
        S_train = np.matmul(df_train[module].values, np.array(coiffiecients).reshape(-1,))
        S_test = np.matmul(df_test[module].values, np.array(coiffiecients).reshape(-1,))

        # Sort the S-scores for training (it's not necessary)
        arg_sort = np.argsort(S_train)
        S_train = S_train[arg_sort]

        # Select the test High and Low risk groups based on the training S-score
        test_HR_indx = S_test > np.percentile(S_train, 66)
        test_LR_indx = S_test < np.percentile(S_train, 33)

        HR_df = df_test.loc[test_HR_indx]
        LR_df = df_test.loc[test_LR_indx]

        results = logrank_test(HR_df['time_to_event'].values, LR_df['time_to_event'].values,
                           HR_df['event'].values, LR_df['event'].values, alpha=.99)
        #print('\n',module)

        stdout = sys.stdout
        sys.stdout = io.StringIO()
        results.print_summary(decimals=5)
        output = sys.stdout.getvalue()
        sys.stdout = stdout


        risk_summary.append('\n'+cancer_+': '+str(module)+'\n'+output)
        # Append the risk analysis to the module's summary
        # with open('../results/CoxPh/{}/Risk_Score_Analysis_summary.txt'.format(m_id), 'a') as f:
        #     f.write(''.join(risk_summary))
