import os
import pandas as pd
import glob
import time

def collect_files(folder):
    '''
    read Genials reports in folder selected, concat the reports in the folder into one dataframe
    :param folder: the folder with Genialis liquid biopsy report to upload, in '.txt' format, get from uploader()
    :return: reports in pandas dataframe
    '''
    all_files = glob.glob(folder + "/*.txt")
    collect = []

    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0, sep="\t")
        df['file_name'] = os.path.basename(filename)
        collect.append(df)

    df = pd.concat(collect, axis=0, ignore_index=True)
    df = df.rename(columns={'Sequence ontology': 'Sequence_ontology', 'Amino Acid Change': 'Amino_Acid_Change',
                            'Allele Frequency': 'Allele_Frequency'})
    df = df.where(pd.notnull(df), None)
    return df

def uploader(replace, conn, folder):
    '''
    FUNCTIONS CHECKS WHETHER ALL FILES FROM FILE_SOURCE TABLE ARE PRESENTD IN MATCHED TABLE (judging criteria is file_name)
    :param replace: IF YOU WISH TO REPLACE CURRRENT FILES PUT replace = True, IF YOU JUST WANT TO ADD SOME NEW FILES WITHOUT REPLACEMENT PUT replace = False
    :param conn: database connector
    :param folder: the folder with Genialis liquid biopsy report to upload, in '.txt' format
    :return: report files to upload
    '''
    # getting unique filenames from the file source table
    datasource_filenames = '''
    select distinct
    file_name
    from FILE_SOURCE
    ;
    '''
    s_fnames = pd.read_sql(datasource_filenames, conn)
    s_fnames = s_fnames.values.tolist()
    s_fnames = [item for sublist in s_fnames for item in sublist]

    files = collect_files(folder)
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
            # print('3. Starting REPLACEMENT in 10 seconds, if you want to CANCEL - stop execution')
            # time.sleep(10)
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
            print('Something went wrong, these files are still not uploaded:', *unmatched_new, sep='\n')
        return unmatched

