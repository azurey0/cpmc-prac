import re
import pandas as pd

##  Amino Acid DICTIONARY, for match in format between Genials reports and COSMIC database
aa_type_dict = {'Ala': 'A', 'Gly': 'G', 'Ile': 'I', 'Leu': 'L', 'Pro': 'P', 'Val': 'V', \
                'Phe': 'F', 'Trp': 'W', 'Tyr': 'Y', 'Asp': 'D', 'Glu': 'E', 'Arg': 'R', 'His': 'H', 'Lys': 'K', \
                'Ser': 'S', 'Thr': 'T', 'Cys': 'C', 'Met': 'M', 'fs': '*', 'Asn': 'N', 'Gln': 'Q'}


def max_aa_mut(mutation_list):
    '''
    find the Amino Acid Change with maximum integers, for example,
    input = [p.Asp502Tyr, p.Asp842Tyr]
    output = p.Asp832Tyr
    :param mutation_list: list of Amino Acid Mutation, corresponding to 'Amino_Acid_Change' values in reports
    :return: Amino Acid Change with maximum integers and it's index in mutation_list
    '''
    values = [re.findall('[0-9]+', x) for x in mutation_list]
    new_list = [item for sublist in values for item in sublist]
    new_list2 = [int(x) for x in new_list]

    return max(new_list2), new_list2.index(max(new_list2))

def checker_filenames(conn):
    '''
    FUNCTIONS CHECKS WHETHER ALL FILES FROM FILE_SOURCE TABLE ARE PRESENTD IN MATCHING_RESULT TABLE
    :param conn: database connection
    :return: unmatched file names
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

    matched_filnames = '''
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


def getting_unmatched_files(conn):
    '''
    execute checker_filenames to get unmatched file names, then get unmatched files' content,
    concat as pandas dataframes
    :param conn: database connection
    :return: unmatched reports concat as dataframe
    '''
    unmatched = checker_filenames(conn)
    unmatched = [item for sublist in unmatched for item in sublist]
    report = pd.DataFrame()
    for i in range(0, len(unmatched)):
        filt = unmatched[i]
        sql_files = f'''
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


def parse_report(conn):
    '''
    Parse report, make report format and tokens usable
    :param conn: database connection
    :return: parsed reports concat as dataframe
    '''
    report = getting_unmatched_files(conn)
    report['Mutation_cosmic'] = ''
    report['Type'] = ''
    report['Score'] = 0

    # RENAME
    for ind, row in report.iterrows():
        value = str(row['Amino_Acid_Change'])
        for key in aa_type_dict.keys():
            value = value.replace(key, aa_type_dict[key])
            report.loc[ind, 'Mutation_cosmic'] = value

    for ind, row in report.iterrows():
        mutation_list = str(row['Mutation_cosmic']).split("|")
        if len(mutation_list) > 1:
            max_score, max_aa_ind = max_aa_mut(mutation_list) # for each row, get max amino acid mutation as its representative
            # print(max_score, max_aa_ind)
            try:
                report.loc[ind, 'Mutation_cosmic'] = mutation_list[max_aa_ind]
            except:
                print('Error: FILE  {}'.format(report.loc[ind, 'file_name']))
                print('Mutation:  {}'.format(report.loc[ind, 'Amino_Acid_Change']))

    return report

def define_right_cosmic(report, conn):
    '''
    An intermediate process to boost running of the program: get mutations from unmatched reports,
    upload to a temp database table for comparison
    :param report: parsed report dataframe
    :param conn: database connection
    :return:
    '''

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


def gene_match(report, cosmic):
    '''
    match gene in reports to COSMIC database by GENE_NAME and MUTATION_AA
    :param report: parsed reports
    :param cosmic: Cosmic dataframe
    :return: reports dataframe with added columns: 'Type' and 'Score', which is matched from COSMIC database
    '''
    #     report = parse_report()
    for ind, row in report.iterrows():
        mut = str(row['Mutation_cosmic'])
        gene_list = str(row['Gene']).split("|")
        for i, gene in enumerate(gene_list):
            # df = cosmic.loc[(cosmic['MUTATION_AA'] == mut) & (cosmic['GENE_NAME'] == gene)]
            df = cosmic[(cosmic['GENE_NAME'] == gene) & (cosmic['MUTATION_AA'] == mut)]
            # print('Searching for: {},{}'.format(gene, mut))   ## here to see what's the current search
            if not df.empty:
                type_1 = df.values[0][3]
                score = df.values[0][4]
                report.loc[ind, 'Type'] = type_1
                report.loc[ind, 'Score'] = score
                # print(gene, mut,type_1,score)

    return report


def matched_table_update(final_result, conn):
    '''
    update MATCHING_RESULT table with processed reports
    :param final_result: processed reports
    :param conn: database connection
    :return:
    '''
    unmatched = checker_filenames(conn)
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

