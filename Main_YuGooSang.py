import time
import os
import multiprocessing as mp
import platform
import random

import Util
import Logic
import LogicPrep
############### st env ################

WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
SYSTEM_NM = platform.system()

if SYSTEM_NM == 'Linux':
    # REAL
    REF_DIR = "/media/backup/ref/GRCm38_p6_mouse/"
else:
    # DEV
    WORK_DIR = "D:/000_WORK/YuGooSang/20201117/WORK_DIR/"
    REF_DIR = "D:/000_WORK/000_reference_path/human/GRCh37_release_101/dna/"

IN = 'input/'
OU = 'output/'

IN_EXCEL = '5&8_input_201127.xlsx'
SHEET_NAME = '8_input_group7_filtered'
############### en env ################


def main():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()

    df = util.read_excel_to_df(WORK_DIR + IN + IN_EXCEL, SHEET_NAME)
    len_df = len(df[df.columns[0]])

    dict_by_chr = {}
    for i in range(len_df):
        dna_cd = df.loc[i][0]
        chr_nm = df.loc[i][1]
        clvg_site = int(df.loc[i][2])
        dna_seq = df.loc[i][3]
        dna_seq_up = df.loc[i][4]
        if chr_nm in dict_by_chr:
            dict_by_chr[chr_nm].append([dna_cd, chr_nm, clvg_site, dna_seq, dna_seq_up])
        else:
            dict_by_chr.update({chr_nm: [[dna_cd, chr_nm, clvg_site, dna_seq, dna_seq_up]]})

    result_list = []
    for ch_key, val_list in dict_by_chr.items():
        p_seq, m_seq = util.read_file_by_biopython(REF_DIR + ch_key + '.fa', 'fasta')

        for val_arr in val_list:
            clvg_site = int(val_arr[2])
            dna_seq_up = val_arr[4]
            p_seq_nr_clvg = p_seq[clvg_site - 18: clvg_site + 5]
            m_seq_nr_clvg = m_seq[clvg_site - 6: clvg_site + 17][::-1]
            if logic.check_strand(p_seq_nr_clvg, m_seq_nr_clvg, dna_seq_up) == '+':
                val_arr.extend(['+', p_seq[clvg_site - 100: clvg_site + 100]])
                result_list.append(val_arr)
            elif logic.check_strand(p_seq_nr_clvg, m_seq_nr_clvg, dna_seq_up) == '-':
                val_arr.extend(['-', m_seq[clvg_site - 101: clvg_site + 99][::-1]])
                result_list.append(val_arr)

    header = ['', 'chr', 'loc', 'seq', 'seq.upper()', 'strand', '100 bp + clvg + 100 bp']
    util.make_excel(WORK_DIR + OU + SHEET_NAME + '_result', header, result_list)


def test():
    """
    + strand
    HEK3-02	chr2	240026760	GGCtCAGACTGAGCACcTGAGAG	GGCTCAGACTGAGCACC TGAGAG
    GGCTCAGACTGAGCACCTGAGAG
    GGCTCAGACTGAGCACCTGAGAG

    - strand
    HEK3-03	chr11	134582414	GGCgCAGACaGAGCACGTGACGA	GGCGCAGACAGAGCACG TGACGA
    GGCGCAGACAGAGCACGTGACGA
    GGCGCAGACAGAGCACGTGACGA

    """

    util = Util.Utils()
    ch_key = 'chr11'
    clvg_site = 134582414
    win0 = 6
    p_seq, m_seq = util.read_file_by_biopython(REF_DIR + ch_key + '.fa', 'fasta')
    print(p_seq[clvg_site - win0 - 12: clvg_site + win0 - 1])
    print(m_seq[clvg_site - win0: clvg_site + win0 + 11])
    print(m_seq[clvg_site - win0: clvg_site + win0 + 11][::-1])


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # test()
    main()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))