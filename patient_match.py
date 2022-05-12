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


def string_cut_1(x):
    '''
    string functions 1
    :param x:
    :return:
    '''
    loc = x.find('-')
    return x[0:loc]


def string_cut_2(x):
    '''
    string functions
    :param x:
    :return:
    '''
    loc = x.find('-ctdna')
    return x[0:loc]


def calculate_sim(left, right):
    '''
    data processing functions for patient matching (the do_magick function)
    :param left:
    :param right:
    :return:
    '''
    left_score = pd.DataFrame(index=left.index, columns=left.index)

    for i in range(left.shape[0]):
        for j in range(left.shape[0]):
            if (j > i):
                left_score.iloc[i, j] = 0
            else:
                left_score.iloc[i, j] = round(cosine_similarity(np.array(left.iloc[i]).reshape(1, -1),
                                                                np.array(left.iloc[j]).reshape(1, -1))[0, 0], 2)

    right_score = pd.DataFrame(index=right.index, columns=right.index)

    for i in range(right.shape[0]):
        for j in range(right.shape[0]):
            if (j > i):
                right_score.iloc[i, j] = 0
            else:
                right_score.iloc[i, j] = round(cosine_similarity(np.array(right.iloc[i]).reshape(1, -1),
                                                                 np.array(right.iloc[j]).reshape(1, -1))[0, 0], 2)
    return left_score, right_score


def calc_sim_insert_index(buckets):
    '''
    data processing functions for patient matching (the do_magick function)
    :param left:
    :param right:
    :return:
    '''
    p_list = []
    p_id = 1
    for bucket in buckets:
        # Define Patient ID
        a = []
        for ind, row in bucket.iterrows():
            b = []
            for rep in range(0, len(row)):
                if ((row[f'{row.index[rep]}'] >= 0.98) and (row.index[rep] != ind)):
                    b.append([row.index[rep]])

            a.append([ind, b])

        # Get the report by patients
        rep_list = []
        for p in range(0, len(a)):
            r_name = a[p][0]
            if len(a[p][1]) == 0:
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
    '''
    data processing functions for patient matching (the do_magick function)
    :param left:
    :param right:
    :return:
    '''
    list_1 = list_1.sort_values()
    max_v = list_1[-2]
    return max_v


def get_scores(buckets):
    '''
    data processing functions for patient matching (the do_magick function)
    :param left:
    :param right:
    :return:
    '''
    collect = []
    for buc in buckets:
        for ind, row in buc.iterrows():
            rep_max = second_largest(row)
            rep_max_ind = row[row == rep_max].index[0]
            collect.append([ind, rep_max, rep_max_ind])
    scores = pd.DataFrame(collect, columns=["file_name", 'max_score', "match"])
    return scores


def get_scores_transposed(buckets):
    '''
    data processing functions for patient matching (the do_magick function)
    :param left:
    :param right:
    :return:
    '''
    collect = []
    for buc in buckets:
        t_buc = buc.T
        for ind, row in t_buc.iterrows():
            rep_max = second_largest(row)
            rep_max_ind = row[row == rep_max].index[0]
            collect.append([ind, rep_max, rep_max_ind])
    scores = pd.DataFrame(collect, columns=["file_name", 'max_score', "match"])
    return scores


def do_magick(data, crc_type):
    '''
    double check matching results.
    :param data:
    :param crc_type: cancer type
    :return: reports with potential matched patients information
    '''
    print(crc_type)
    # Data pre processing
    data_cancer = data[(data['cancer_type'] == crc_type)]
    data_cancer_pivot = pd.pivot_table(data_cancer, values='Allele_Frequency',
                                       index='file_name_sc', columns='Gene',
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
    condR = tempR.iloc[3] < 0.8
    data_cancer_pivot_R = data_cancer_pivot_R.loc[:, condR]

    tempL = pd.DataFrame(data_cancer_pivot_L.describe())
    condL = tempL.iloc[3] < 0.8
    data_cancer_pivot_L = data_cancer_pivot_L.loc[:, condL]

    print(data_cancer_pivot_L.shape)
    print(data_cancer_pivot_R.shape)

    # 1st Round

    a, b = calculate_sim(data_cancer_pivot_L, data_cancer_pivot_R)
    buckets = [a, b]
    p_list = calc_sim_insert_index(buckets)

    df_L = p_list[0]
    l_i = df_L[df_L['Matched'] == False]
    l_m = df_L[df_L['Matched'] == True]

    df_R = p_list[1]
    r_i = df_R[df_R['Matched'] == False]
    r_m = df_R[df_R['Matched'] == True]

    ll_list = [l_m, r_i]
    ll_df = pd.concat(ll_list, axis=0, ignore_index=False).reset_index(drop=True)

    rr_list = [r_m, l_i]
    rr_df = pd.concat(rr_list, axis=0, ignore_index=False).reset_index(drop=True)

    data_cancer_pivot_LL = data_cancer_pivot[data_cancer_pivot.index.isin(ll_df['file_name'])]

    data_cancer_pivot_RR = data_cancer_pivot[data_cancer_pivot.index.isin(rr_df['file_name'])]

    # 2nd Round

    c, d = calculate_sim(data_cancer_pivot_LL, data_cancer_pivot_RR)
    buckets_1 = [c, d]
    p_list_1 = calc_sim_insert_index(buckets_1)

    df_LL = p_list_1[0]
    ll_i = df_LL[df_LL['Matched'] == False]
    ll_m = df_LL[df_LL['Matched'] == True]

    df_RR = p_list_1[1]
    rr_i = df_RR[df_RR['Matched'] == False]
    rr_m = df_RR[df_RR['Matched'] == True]

    matched_list = [ll_m, rr_m]
    matched_df = pd.concat(matched_list, axis=0, ignore_index=False).reset_index(drop=True)

    individual_list = [rr_i, ll_i]
    individual_df = pd.concat(individual_list, axis=0, ignore_index=False).reset_index(drop=True)

    data_cancer_pivot_matched = data_cancer_pivot[data_cancer_pivot.index.isin(matched_df['file_name'])]
    data_cancer_pivot_individual = data_cancer_pivot[data_cancer_pivot.index.isin(individual_df['file_name'])]

    # Last Round

    e, f = calculate_sim(data_cancer_pivot_matched, data_cancer_pivot_individual)
    buckets_final = [e, f]

    p_list_final = calc_sim_insert_index(buckets_final)
    final_patients = pd.concat(p_list_final, axis=0, ignore_index=False).reset_index(drop=True)
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
    final_patients['Best_round'] = final_patients[['score_r1', 'score_r2', 'score_r3']].idxmax(axis=1)

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
        if len(a) > 0:
            final_patients.loc[ind, 'Potential_patient'] = a.Patient_ID.values[0]
        else:
            final_patients.loc[ind, 'Potential_patient'] = 'None'

    df_out = final_patients[['Cancer_Type', 'file_name', 'Matched',
                             'Patient_ID', 'Best_Score', 'Potential_patient']].reset_index(drop=True)

    df_out = df_out.rename(columns={"file_name": "file_name_sc"})

    return df_out


def patient_table_update(patients_to_update, conn):
    '''
    update PATIENT table
    :param patients_to_update: reports with potential matched patient information, like the 'df_out' result above
    :param conn: database connection
    :return:
    '''
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

