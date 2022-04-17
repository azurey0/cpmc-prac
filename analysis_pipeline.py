import file_uploader
import gene_match
import patient_match
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

def folder_check(folder):
    error_msg = None
    # check if folder exists
    if os.path.isdir(folder)==False:
        error_msg = 'Folder does not exists.'
    else:
        # check if folder contains txt files
        all_files = glob.glob(folder + "/*.txt")
        if not all_files:
            error_msg = 'No .txt file exists in the folder.'
        else:
            # check if file contains important columns
            for file in all_files:
                df = pd.read_csv(file, index_col=None, header=0, sep="\t")
                col_check = set(df.columns).issuperset(['Sequence ontology','Amino Acid Change','Allele Frequency'])
                if col_check == False:
                    error_msg = 'File does not contain either Sequence ontology, Amino Acid Change or Allele Frequency.'
    return error_msg



def report_analysis(folder, replace=True):
    '''
    gui touchpoint
    :param foler: the folder with liquid biopsy reports to analyze.
    :param replace: whether it is to replace(update) existing reports or not
    :return: database updates.
    '''
    conn = pymysql.connect(
        host='cpmc.ccwhfsqzn9mq.us-west-1.rds.amazonaws.com',
        port=int(3306),
        user="admin",
        passwd="UCDavis2021",
        db="cpmcpracticum",
        charset='utf8mb4')

    print('Succesfully connected!')

    # 1

    # Execute files upload

    # IF YOU WISH TO REPLACE CURRRENT FILES PUT replace = True

    # IF YOU JUST WANT TO ADD SOME NEW FILES WITHOUT REPLACEMENT PUT replace = False

    #folder_path = 'C:\\Users\\92350\\PycharmProjects\\gui\\DATA_CPMC'
    #foler_path = folder
    start = time.time()
    unmatched = file_uploader.uploader(replace=replace, conn=conn, folder=folder)
    end = time.time()
    print('update took {} sec'.format(round(end-start,2)))


    # 2

    # Execute GENE MATCHING

    # Creating the base version of a report
    report = gene_match.parse_report(conn)

    start = time.time()

    gene_match.define_right_cosmic(report, conn)
    # COSMIC DF
    get_cosmic='''
    select 
    *
    from COSMIC c 
    where 1=1
    and c.MUTATION_AA in(select
                         mtm.Mutation_cosmic
                         from mut_to_match mtm)
    ;
    '''

    cosmic = pd.read_sql(get_cosmic, conn)

    end = time.time()
    print(round(end-start,2))

    # GENERATING FINAL RESULT
    report_result = gene_match.gene_match(report, cosmic)

    ## write the result to appropriate format
    #today = date.today()
    #d = today.strftime("%b-%d-%Y")
    final_result = pd.DataFrame(report_result,columns=['file_name','Gene', 'Amino_Acid_Change',
                                                       'Allele_Frequency','Mutation_cosmic','Type','Score'])
    final_result = final_result.where(pd.notnull(final_result), None)
    # final_result.to_csv('matching_result_{}.csv'.format(d))

    # DATABASE UPDATE
    gene_match.matched_table_update(final_result, conn)

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

    data['cancer_type'] = data.apply(lambda x: patient_match.string_cut_1(x['file_name']), axis = 1)
    data['file_name_sc'] = data.apply(lambda x: patient_match.string_cut_2(x['file_name']), axis = 1)

    df_out = data[['file_name', 'file_name_sc']].drop_duplicates().reset_index(drop=True)


    cncr_types = data['cancer_type'].unique()

    # Set Insturctions
    new = 0

    collect = []
    if new == 0:
        for crc in cncr_types:
            try:
                crc_df = patient_match.do_magick(data, crc)
                collect.append(crc_df)
            except:
                pass
        all_reports = pd.concat(collect, axis=0, ignore_index = False).reset_index(drop = True)

    df_out = df_out.merge(all_reports, how='inner', on='file_name_sc')
    df_out['Matched'] = df_out['Matched'].astype(int)
    df_out = df_out.where(pd.notnull(df_out), None)

    # Update Database
    patient_match.patient_table_update(df_out, conn)


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
    t3.file_name_sc,
    t3.Matched,
    t3.Patient_ID,
    t3.Best_Score,
    t3.Potential_patient,
    pm.num_match
    FROM t2
    JOIN t1
    left join t3 on t2.file_name = t3.file_name
    left join (select distinct
                p.Patient_ID,
                count(p.file_name) as num_match
                from PATIENTS p
                group by 1
                ) pm on t3.Patient_ID = pm.Patient_ID
    WHERE 1=1
    and t1.file_name = t2.file_name
    ;
    '''

    cursor.execute(sql_DB_table)
    conn.commit()

### TESTING
#folder_path = 'C:\\Users\\92350\\PycharmProjects\\gui\\DATA_CPMC'
#report_analysis(folder_path, 1)


