### PACKAGES##
import os
import pandas as pd
import pymysql
import glob
import time
from datetime import date
import re
from tqdm import tqdm
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import cosine_similarity,cosine_distances
import numpy as np
import seaborn as sns


conn = pymysql.connect(
    host='cpmc.ccwhfsqzn9mq.us-west-1.rds.amazonaws.com',
    port=int(3306),
    user="admin",
    passwd="UCDavis2021",
    db="cpmcpracticum",
    charset='utf8mb4')

print('Succesfully connected!')


### FUNCTIONS ###


############### Files Upload
def collect_files(file_path):
    path = os.path.join(os.getcwd(),'DATA_CPMC')
    all_files = glob.glob(path + "/*.txt")
    collect = []

    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0, sep="\t")
        df['file_name'] = os.path.basename(filename)
        collect.append(df)

    df = pd.concat(collect, axis=0, ignore_index=True)
    df = df.rename(columns = {'Sequence ontology':'Sequence_ontology', 'Amino Acid Change':'Amino_Acid_Change', 'Allele Frequency':'Allele_Frequency'})
    df = df.where(pd.notnull(df), None)
    print(df[:5])
    return df

# FUNCTIONS CHECKS WHETHER ALL FILES FROM DS TABLE ARE PRESENTD IN MATCHED TABLE
def uploader(replace):
    # getting unique filenames from the file source table
    datasource_filenames='''
    select distinct
    file_name
    from FILE_SOURCE
    ;
    '''
    s_fnames = pd.read_sql(datasource_filenames, conn)
    s_fnames = s_fnames.values.tolist()
    s_fnames = [item for sublist in s_fnames for item in sublist]
    
    files = collect_files()
    files_to_upload = files['file_name'].unique()
    files_to_upload = files_to_upload.tolist()
    
    unmatched = list()
    for f in files_to_upload:
        if f not in s_fnames:
            unmatched.append(f)
    if len(unmatched) == 0:
        print('1. Check is complete: NO NEW FILES FOUND')
        if replace == True:
            print('2. updating database...')
            # INSERT all reports with replacement
            print('3. Starting REPLACEMENT in 10 seconds, if you want to CANCEL - stop execution')
            time.sleep(10)
            print('4. executing...')
            cursor = conn.cursor()
            for i in range(0, len(files_to_upload)):
                filt = files_to_upload[i]
                deleting_stuff = f'''delete
                                        from FILE_SOURCE
                                        where 1=1
                                        and file_name in('{filt}')'''
                cursor.execute(deleting_stuff)
                conn.commit()
            
            # INSERT VALUES of the result
            values = list(map(tuple, files.values))
            cursor.executemany('''INSERT INTO FILE_SOURCE(
            Chromosome,
            Position,
            ID,
            Impact,
            Sequence_ontology,
            Amino_Acid_Change,
            Gene,
            Allele_Frequency,
            file_name
            ) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)
            ''', values)
            conn.commit()
            print('Database Updated')
            return unmatched
        else:
            print('Database not updated')
    else:
        print(f'1. FOUND {len(unmatched)} NEW FILES')
        print('2. updating database...')
        if replace == True:
            # INSERT all reports with replacement
            print('3. Starting REPLACEMENT in 10 seconds, if you want to CANCEL - stop execution')
            time.sleep(10)
            print('4. executing...')
            cursor = conn.cursor()
            for i in range(0, len(files_to_upload)):
                filt = files_to_upload[i]
                deleting_stuff = f'''delete
                                        from FILE_SOURCE
                                        where 1=1
                                        and file_name in('{filt}')'''
                cursor.execute(deleting_stuff)
                conn.commit()
            # INSERT VALUES of the result
            values = list(map(tuple, files.values))
            cursor.executemany('''INSERT INTO FILE_SOURCE(
            Chromosome,
            Position,
            ID,
            Impact,
            Sequence_ontology,
            Amino_Acid_Change,
            Gene,
            Allele_Frequency,
            file_name
            ) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)
            ''', values)
            conn.commit()
        else:
            # INSERT only reports that are not currently present without replacement
            files_no_replacement = files[files['file_name'].isin(unmatched)]
            cursor = conn.cursor()
            values = list(map(tuple, files_no_replacement.values))
            cursor.executemany('''INSERT INTO FILE_SOURCE(
            Chromosome,
            Position,
            ID,
            Impact,
            Sequence_ontology,
            Amino_Acid_Change,
            Gene,
            Allele_Frequency,
            file_name
            ) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)
            ''', values)
            conn.commit()
            
        s_fnames_new = pd.read_sql(datasource_filenames, conn)
        s_fnames_new = s_fnames_new.values.tolist()
        s_fnames_new = [item for sublist in s_fnames_new for item in sublist]
        unmatched_new = list()
        for f in files_to_upload:
            if f not in s_fnames_new:
                unmatched_new.append(f)
        if len(unmatched_new) == 0:
            print('Database Sucessfully Updated!')
        else:
            print('Something went wrong, these files are still not uploaded:', *unmatched_new, sep = '\n')
        return unmatched

    
    
############### Matching Genes

## DICTIONARY
aa_type_dict={'Ala':'A', 'Gly':'G', 'Ile':'I', 'Leu':'L', 'Pro':'P', 'Val':'V',\
        'Phe':'F', 'Trp':'W', 'Tyr':'Y', 'Asp':'D', 'Glu':'E', 'Arg':'R', 'His':'H', 'Lys':'K',\
        'Ser':'S', 'Thr':'T', 'Cys':'C', 'Met':'M', 'fs':'*', 'Asn' : 'N', 'Gln':'Q'}

def max_aa_mut(mutation_list):
    values = [re.findall('[0-9]+', x) for x in mutation_list]
    new_list = [item for sublist in values for item in sublist]
    new_list2 = [int(x) for x in new_list]
    
    return max(new_list2), new_list2.index(max(new_list2))

# FUNCTIONS CHECKS WHETHER ALL FILES FROM DS TABLE ARE PRESENTD IN MATCHED TABLE
def checker_filenames():
    # getting unique filenames from the file source table
    datasource_filenames='''
    select distinct
    file_name
    from FILE_SOURCE
    ;
    '''
    s_fnames = pd.read_sql(datasource_filenames, conn)
    s_fnames = s_fnames.values.tolist()
    
    matched_filnames='''
    select distinct
    file_name
    from MATCHING_RESULT
    ;
    '''
    m_fnames = pd.read_sql(matched_filnames, conn)
    m_fnames = m_fnames.values.tolist()
    
    unmatched = list()
    for f in s_fnames:
        if f not in m_fnames:
            unmatched.append(f)
    if len(unmatched) == 0:
        print('Check is complete: all current files are already matched!')
        return unmatched
    else:
        print('FOUND UNMATCHED FILES, processing...')
        return unmatched

def getting_unmatched_files():
    unmatched = checker_filenames()
    unmatched = [item for sublist in unmatched for item in sublist]
    report = pd.DataFrame()
    for i in range(0, len(unmatched)):
        filt = unmatched[i]
        sql_files=f'''
        select 
        *
        from FILE_SOURCE
        where 1=1
        and file_name in('{filt}')
        ;
        '''
        file = pd.read_sql(sql_files, conn)
        report = pd.concat([report, file], axis=0).reset_index(drop=True)
    return report

def parse_report():
    report = getting_unmatched_files()
    report['Mutation_cosmic'] = ''
    report['Type'] = ''
    report['Score'] = 0
    
    # RENAME
    for ind, row in tqdm(report.iterrows(), total=report.shape[0]):
        value = str(row['Amino_Acid_Change'])
        for key in aa_type_dict.keys():
            value = value.replace(key, aa_type_dict[key])
            report.loc[ind, 'Mutation_cosmic'] = value
    
    for ind, row in tqdm(report.iterrows(), total=report.shape[0]):    
        mutation_list = str(row['Mutation_cosmic']).split("|")
        if len(mutation_list) > 1:
            max_score, max_aa_ind = max_aa_mut(mutation_list)
            #print(max_score, max_aa_ind)
            try:
                report.loc[ind, 'Mutation_cosmic'] = mutation_list[max_aa_ind]
            except:
                print('Error: FILE  {}'.format(report.loc[ind, 'file_name']))
                print('Mutation:  {}'.format(report.loc[ind, 'Amino_Acid_Change']))
            
    
    return report

def define_right_cosmic(report):
#     report = parse_report()
    mutations = pd.DataFrame(report['Mutation_cosmic'].unique(), columns = ['Mutation_cosmic'])
    cursor = conn.cursor()
    truncate = '''
    truncate table mut_to_match
    ;
    '''
    cursor.execute(truncate)
    conn.commit()
    values = list(map(tuple, mutations.values))
    cursor.executemany('''INSERT INTO mut_to_match(
                            Mutation_cosmic
                            ) VALUES (%s)
                            ;
                        ''', values)
    conn.commit()
    return
    
def gene_match(report):
#     report = parse_report()
    for ind, row in tqdm(report.iterrows(), total=report.shape[0]):  
        mut = str(row['Mutation_cosmic'])
        gene_list = str(row['Gene']).split("|")
        for i,gene in enumerate(gene_list):      
            # df = cosmic.loc[(cosmic['MUTATION_AA'] == mut) & (cosmic['GENE_NAME'] == gene)]
            df = cosmic[(cosmic['GENE_NAME']== gene) & (cosmic['MUTATION_AA']== mut)]
            #print('Searching for: {},{}'.format(gene, mut))   ## here to see what's the current search
            if not df.empty:
                type_1 = df.values[0][3]
                score = df.values[0][4]
                report.loc[ind, 'Type'] = type_1
                report.loc[ind, 'Score'] = score
                #print(gene, mut,type_1,score)
        
    return report

def matched_table_update(final_result):
    unmatched = checker_filenames()
    # FIRST DELETE IF THEY WERE PRESENT IN THE MATCHED TABLE
    cursor = conn.cursor()
    for i in range(0, len(unmatched)):
        filt = unmatched[i][0]
        deleting_stuff = f'''delete
                                from MATCHING_RESULT
                                where 1=1
                                and file_name in('{filt}')'''
        cursor.execute(deleting_stuff)
        
    conn.commit()
    # INSERT VALUES of the result
    values = list(map(tuple, final_result.values))
    cursor.executemany('''INSERT INTO MATCHING_RESULT(
    file_name,
    Gene,
    Amino_Acid_Change,
    Allele_Frequency,
    Mutation_cosmic,
    Type,
    Score
    ) VALUES (%s,%s,%s,%s,%s,%s,%s)
    ''', values)
    conn.commit()
    return



################ PATIENT MATCHING

def string_cut_1(x):
    loc = x.find('-')
    return x[0:loc]


def string_cut_2(x):
    loc = x.find('-ctdna')
    return x[0:loc]


def calculate_sim(left, right):
    
    left_score = pd.DataFrame(index=left.index, columns = left.index)

    for i in tqdm(range(left.shape[0])):
        for j in range(left.shape[0]):
            if (j>i):
                left_score.iloc[i,j] = 0
            else:
                left_score.iloc[i,j]=round(cosine_similarity(np.array(left.iloc[i]).reshape(1,-1),
                                                             np.array(left.iloc[j]).reshape(1,-1))[0,0],2)

    right_score = pd.DataFrame(index=right.index, columns = right.index)

    for i in tqdm(range(right.shape[0])):
        for j in range(right.shape[0]):
            if (j>i):
                right_score.iloc[i,j] = 0
            else:
                right_score.iloc[i,j]=round(cosine_similarity(np.array(right.iloc[i]).reshape(1,-1),
                                                             np.array(right.iloc[j]).reshape(1,-1))[0,0],2)
    return left_score, right_score


def calc_sim_insert_index(buckets):
    p_list = []
    p_id = 1
    for bucket in buckets:
        # Define Patient ID
        a = []
        for ind, row in bucket.iterrows():
            b = []
            for rep in range(0,len(row)):
                if ((row[f'{row.index[rep]}'] >=0.98) and (row.index[rep] != ind)):
                    b.append([row.index[rep]])

                        
            a.append([ind,b])

        # Get the report by patients
        rep_list = []
        for p in range(0,len(a)):
            r_name = a[p][0]
            if len(a[p][1])==0:
                patient_id = p_id
                p_id += 1
            else:
                key = a[p][1][0][0]
                val = [elt[1] for elt in rep_list if elt[0] == key][0]
                patient_id = val

            rep_list.append([r_name, patient_id])
        df = pd.DataFrame(rep_list, columns=["file_name", "patient"])
        df['Matched'] = df['patient'].duplicated(keep=False)
        p_list.append(df)
    return p_list


def second_largest(list_1):
    list_1 = list_1.sort_values()
    max_v = list_1[-2]
    return max_v

def get_scores(buckets):
    collect = []
    for buc in buckets:
        for ind, row in buc.iterrows():  
            rep_max = second_largest(row)
            rep_max_ind = row[row == rep_max].index[0]
            collect.append([ind, rep_max, rep_max_ind])
    scores = pd.DataFrame(collect, columns=["file_name",'max_score', "match"])
    return scores

def get_scores_transposed(buckets):
    collect = []
    for buc in buckets:
        t_buc = buc.T
        for ind, row in t_buc.iterrows():  
            rep_max = second_largest(row)
            rep_max_ind = row[row == rep_max].index[0]
            collect.append([ind, rep_max, rep_max_ind])
    scores = pd.DataFrame(collect, columns=["file_name",'max_score', "match"])
    return scores


def do_magick(crc_type):
    print(crc_type)
    # Data pre processing
    data_cancer = data[(data['cancer_type'] == crc_type)]
    data_cancer_pivot = pd.pivot_table(data_cancer, values='Allele_Frequency',
                                      index = 'file_name_sc',columns = 'Gene',
                                      aggfunc=max, fill_value=0)
    
    
    # Split into buckets
    X = data_cancer_pivot
    km = KMeans(n_clusters=2, init='random')
    y_km = km.fit_predict(X)
    data_cancer_pivot['Bucket'] = y_km
    data_cancer_pivot_L = data_cancer_pivot.loc[data_cancer_pivot['Bucket'] == 1]
    data_cancer_pivot_L = data_cancer_pivot_L.loc[:, data_cancer_pivot_L.columns != 'Bucket']
    data_cancer_pivot_R = data_cancer_pivot.loc[data_cancer_pivot['Bucket'] == 0]
    data_cancer_pivot_R = data_cancer_pivot_R.loc[:, data_cancer_pivot_R.columns != 'Bucket']
    
    # Get rid of some mutations
    tempR = pd.DataFrame(data_cancer_pivot_R.describe())
    condR = tempR.iloc[3]<0.8
    data_cancer_pivot_R = data_cancer_pivot_R.loc[:, condR]
    
    tempL = pd.DataFrame(data_cancer_pivot_L.describe())
    condL = tempL.iloc[3]<0.8
    data_cancer_pivot_L = data_cancer_pivot_L.loc[:, condL]
    
    print(data_cancer_pivot_L.shape)
    print(data_cancer_pivot_R.shape)
    
    # 1st Round
    
    a,b = calculate_sim(data_cancer_pivot_L, data_cancer_pivot_R)
    buckets = [a,b]
    p_list = calc_sim_insert_index(buckets)
    
    df_L = p_list[0]
    l_i = df_L[df_L['Matched'] == False]
    l_m = df_L[df_L['Matched'] == True]

    df_R = p_list[1]
    r_i = df_R[df_R['Matched'] == False]
    r_m = df_R[df_R['Matched'] == True]

    ll_list = [l_m, r_i]
    ll_df = pd.concat(ll_list, axis=0, ignore_index = False).reset_index(drop = True)

    rr_list = [r_m, l_i]
    rr_df = pd.concat(rr_list, axis=0, ignore_index = False).reset_index(drop = True)

    data_cancer_pivot_LL = data_cancer_pivot[data_cancer_pivot.index.isin(ll_df['file_name'])]

    data_cancer_pivot_RR = data_cancer_pivot[data_cancer_pivot.index.isin(rr_df['file_name'])]
    
    # 2nd Round
    
    c,d = calculate_sim(data_cancer_pivot_LL, data_cancer_pivot_RR)
    buckets_1 = [c,d]
    p_list_1 = calc_sim_insert_index(buckets_1)
    
    df_LL = p_list_1[0]
    ll_i = df_LL[df_LL['Matched'] == False]
    ll_m = df_LL[df_LL['Matched'] == True]

    df_RR = p_list_1[1]
    rr_i = df_RR[df_RR['Matched'] == False]
    rr_m = df_RR[df_RR['Matched'] == True]

    matched_list = [ll_m, rr_m]
    matched_df = pd.concat(matched_list, axis=0, ignore_index = False).reset_index(drop = True)

    individual_list = [rr_i, ll_i]
    individual_df = pd.concat(individual_list, axis=0, ignore_index = False).reset_index(drop = True)

    data_cancer_pivot_matched = data_cancer_pivot[data_cancer_pivot.index.isin(matched_df['file_name'])]
    data_cancer_pivot_individual = data_cancer_pivot[data_cancer_pivot.index.isin(individual_df['file_name'])]
    
    # Last Round

    e,f = calculate_sim(data_cancer_pivot_matched, data_cancer_pivot_individual)
    buckets_final = [e,f]
    
    p_list_final = calc_sim_insert_index(buckets_final)
    final_patients = pd.concat(p_list_final, axis=0, ignore_index = False).reset_index(drop = True)
    final_patients['Cancer_Type'] = crc_type
    final_patients['Patient_ID'] = final_patients['Cancer_Type'] + '-' + final_patients['patient'].astype(str)
    
    # Getting Scores and Flags
    
    # Round 1
    r1_scores_n = get_scores(buckets)
    r1_scores_t = get_scores_transposed(buckets)
    r1_scores = r1_scores_n.merge(r1_scores_t, how='inner', on='file_name')
    r1_scores['score_r1'] = 0
    r1_scores['best_r1'] = ''
    for ind, row in r1_scores.iterrows():
        if row['max_score_x'] >= row['max_score_y']:
            r1_scores.loc[ind, 'score_r1'] = row['max_score_x']
            r1_scores.loc[ind, 'best_r1'] = row['match_x']
        else:
            r1_scores.loc[ind, 'score_r1'] = row['max_score_y']
            r1_scores.loc[ind, 'best_r1'] = row['match_y']
    r1_scores = r1_scores[['file_name', 'score_r1', 'best_r1']]
    
    # Round 2
    r2_scores_n = get_scores(buckets_1)
    r2_scores_t = get_scores_transposed(buckets_1)
    r2_scores = r2_scores_n.merge(r2_scores_t, how='inner', on='file_name')
    r2_scores['score_r2'] = 0
    r2_scores['best_r2'] = ''
    for ind, row in r2_scores.iterrows():
        if row['max_score_x'] >= row['max_score_y']:
            r2_scores.loc[ind, 'score_r2'] = row['max_score_x']
            r2_scores.loc[ind, 'best_r2'] = row['match_x']
        else:
            r2_scores.loc[ind, 'score_r2'] = row['max_score_y']
            r2_scores.loc[ind, 'best_r2'] = row['match_y']
    r2_scores = r2_scores[['file_name', 'score_r2', 'best_r2']]
    
    # Round 3
    r3_scores_n = get_scores(buckets_final)
    r3_scores_t = get_scores_transposed(buckets_final)
    r3_scores = r3_scores_n.merge(r3_scores_t, how='inner', on='file_name')
    r3_scores['score_r3'] = 0
    r3_scores['best_r3'] = ''
    for ind, row in r3_scores.iterrows():
        if row['max_score_x'] >= row['max_score_y']:
            r3_scores.loc[ind, 'score_r3'] = row['max_score_x']
            r3_scores.loc[ind, 'best_r3'] = row['match_x']
        else:
            r3_scores.loc[ind, 'score_r3'] = row['max_score_y']
            r3_scores.loc[ind, 'best_r3'] = row['match_y']
    r3_scores = r3_scores[['file_name', 'score_r3', 'best_r3']]
    
    # Getting all together
    final_patients = final_patients.merge(r1_scores, how='inner', on='file_name')
    final_patients = final_patients.merge(r2_scores, how='inner', on='file_name')
    final_patients = final_patients.merge(r3_scores, how='inner', on='file_name')
    final_patients['Best_round'] = final_patients[['score_r1', 'score_r2', 'score_r3']].idxmax(axis = 1)
    
    final_patients = final_patients.sort_values(by=['Matched', 'patient'])
    
    final_patients['Best_Score'] = 0
    final_patients['Best_match'] = ''
    final_patients['Potential_patient'] = ''
    
    for ind, row in final_patients.iterrows():
        b_round = row['Best_round']
        final_patients.loc[ind, 'Best_Score'] = row[f'{b_round}']
    
    for ind, row in final_patients.iterrows():
        if row['Best_Score'] >= 0.9:
            if row['Best_round'] == 'score_r1':
                final_patients.loc[ind, 'Best_match'] = row['best_r1']
            elif row['Best_round'] == 'score_r2':
                final_patients.loc[ind, 'Best_match'] = row['best_r2']
            else:
                final_patients.loc[ind, 'Best_match'] = row['best_r3']
        else:
            final_patients.loc[ind, 'Best_match'] = 'None'
    
    
    for ind, row in final_patients.iterrows():
        a = final_patients[final_patients['file_name'] == row['Best_match']]
        if len(a)>0:
            final_patients.loc[ind, 'Potential_patient'] = a.Patient_ID.values[0]
        else:
            final_patients.loc[ind, 'Potential_patient'] = 'None'
    
    
    df_out = final_patients[['Cancer_Type', 'file_name', 'Matched',
                             'Patient_ID','Best_Score', 'Potential_patient']].reset_index(drop = True)
    
    df_out = df_out.rename(columns={"file_name": "file_name_sc"})
    
    return df_out

def patient_table_update(patients_to_update):
    cursor = conn.cursor()
# Clear the table
    sql1 = """
    truncate table PATIENTS;
    """
    cursor.execute(sql1)
    conn.commit()
    
# Insert new values
    values = list(map(tuple, patients_to_update.values))
    cursor.executemany('''INSERT INTO PATIENTS(
    file_name,
    file_name_sc,
    Cancer_Type,
    Matched,
    Patient_ID,
    Best_Score,
    Potential_patient
    ) VALUES (%s,%s,%s,%s,%s,%s,%s)
    ''', values)
    conn.commit()
    return


################# EXECUTION #################

# 1

# Execute files upload

# IF YOU WISH TO REPLACE CURRRENT FILES PUT replace = True

# IF YOU JUST WANT TO ADD SOME NEW FILES WITHOUT REPLACEMENT PUT replace = False

start = time.time()
unmatched = uploader(replace = False)
end = time.time()
print('update took {} sec'.format(round(end-start,2)))


# 2

# Execute GENE MATCHING

# Creating the base version of a report
report = parse_report()

start = time.time()

define_right_cosmic(report)
# COSMIC DF
sql='''
select 
*
from COSMIC c 
where 1=1
and c.MUTATION_AA in(select
                     mtm.Mutation_cosmic
                     from mut_to_match mtm)
;
'''

cosmic = pd.read_sql(sql, conn)
cosmic[:5]
cosmic.shape

end = time.time()
print(round(end-start,2))

# GENERATING FINAL RESULT
report_result = gene_match(report)

## write the result to appropriate format and to .csv file
today = date.today()
d = today.strftime("%b-%d-%Y")
final_result = pd.DataFrame(report_result,columns=['file_name','Gene', 'Amino_Acid_Change',
                                                   'Allele_Frequency','Mutation_cosmic','Type','Score'])
final_result = final_result.where(pd.notnull(final_result), None)
# final_result.to_csv('matching_result_{}.csv'.format(d))

# DATABASE UPDATE
matched_table_update(final_result)


# 3

# PATIENT MATCHING

sql_matched = '''
select
*
from MATCHING_RESULT
where 1=1
'''

data = pd.read_sql(sql_matched,conn)
# or
# data = final_result

data['cancer_type'] = data.apply(lambda x: string_cut_1(x['file_name']), axis = 1)
data['file_name_sc'] = data.apply(lambda x: string_cut_2(x['file_name']), axis = 1)

df_out = data[['file_name', 'file_name_sc']].drop_duplicates().reset_index(drop=True)


cncr_types = data['cancer_type'].unique()

# Set Insturctions
new = 0

collect = []
if new == 0:
    for crc in cncr_types:
        try:
            crc_df = do_magick(crc)
            collect.append(crc_df)
        except:
            pass
    all_reports = pd.concat(collect, axis=0, ignore_index = False).reset_index(drop = True)

df_out = df_out.merge(all_reports, how='inner', on='file_name_sc')
df_out['Matched'] = df_out['Matched'].astype(int) 
df_out = df_out.where(pd.notnull(df_out), None)

# Update Database
patient_table_update(df_out)


# 4

# DASHBOARD TABLE CREATION

cursor = conn.cursor()
# Clear the table
sql1 = """
drop table if exists tableau_data_2
"""
cursor.execute(sql1)
conn.commit()


sql_DB_table = '''
CREATE TABLE tableau_data_2 AS

with
-- Query for number of distinct results per file
t1 as
(SELECT
 file_name,
 count(DISTINCT GENE) AS num_genes
 FROM MATCHING_RESULT
 group by file_name
 order by count(GENE)),

t2 as
(SELECT *,

-- Pull out cancer type
-- SUBSTRING_INDEX(file_name, "-", 1) AS cancer_type,

-- Germline identification.
CASE
    WHEN Allele_Frequency >= 0.4 THEN 'Germline'
    ELSE ''
END AS Germline_ID,

-- Query for af threshold identification
CASE
    WHEN Allele_Frequency <= 0.01 THEN 'Low_AF'
    ELSE 'High_AF'
END AS AF_Threshold,

-- Query for number of results per file
count(GENE) OVER (PARTITION BY file_name) AS num_results
FROM MATCHING_RESULT),

-- Query for patients & cancer type
t3 as
(SELECT
file_name,
file_name_sc,
Cancer_Type,
Matched,
Patient_ID,
Best_Score,
Potential_patient
FROM PATIENTS)

SELECT
t3.Cancer_Type as cancer_type,
t2.*,
t1.num_genes,
-- Calculate ratio
t1.num_genes/t2.num_results AS distinct_ratio,
-- Flag for bad files, update 32 & 0.25 thresholds if needed
CASE
    WHEN num_results <= 32 THEN 'Low Count'
    WHEN t1.num_genes/t2.num_results <= 0.25 THEN 'Low Ratio'
    ELSE ''
END AS Bad_File_Flag,
t3.Matched,
t3.Patient_ID,
t3.Best_Score,
t3.Potential_patient
FROM t2
JOIN t1
left join t3 on t2.file_name = t3.file_name
WHERE 1=1
and t1.file_name = t2.file_name
;
'''

cursor.execute(sql_DB_table)
conn.commit()