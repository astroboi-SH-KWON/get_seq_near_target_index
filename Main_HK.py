import time
import os
from Bio import SeqIO
import multiprocessing as mp
import numpy as np
import platform

import Util
import Logic
import LogicPrep
############### start to set env ################
WORK_DIR = os.getcwd() + "/"
SYSTEM_NM = platform.system()
if SYSTEM_NM == 'Linux':
    # REAL
    REF_DIR = "/media/hkim/Pipeline/project_by_astroboi_2200/hg38/"
else:
    # DEV
    REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"

PROJECT_NAME = WORK_DIR.split("/")[-2]
CDS_INFO = "hg38_refFlat_full.txt"
MUT_INFO = "20200909_ClinVar_hg38.txt"

PAM = 'NGG'
WIN_PAM = [10, 10]
SEQ_WIN_SIZE = [30, 30]


TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
############### end setting env #################

def multi_processing():
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    # CHROM	POS	ID	REF	ALT	mut_length	CLNVC	CLNSIG
    # POS - 1 = index of .fa sequence
    # [['1', '930188', '846933', 'G', 'A', '1', 'substitution', 'Uncertain_significance'],...]
    mut_list = []
    if SYSTEM_NM == 'Linux':
        mut_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t"))
    else:
        mut_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t")[:300])

    splited_mut_list = np.array_split(mut_list, MULTI_CNT)

    print("platform.system() : ", SYSTEM_NM)
    print("total cpu_count : ", str(TOTAL_CPU))
    print("will use : ", str(MULTI_CNT))
    pool = mp.Pool(processes=MULTI_CNT)

    pool_list = pool.map(get_PAM_within_N_bp_of_POS, splited_mut_list)
    result_list = logic_prep.merge_multi_list(pool_list)

    header = ['CHROM', 'PAM', str(SEQ_WIN_SIZE[0]) + ' + PAM + ' + str(SEQ_WIN_SIZE[1]), 'PAM_POS', 'STRAND']
    util.make_excel(WORK_DIR + "output/ClinVar_hg38_result", header, result_list)

def get_PAM_within_N_bp_of_POS(mut_list):
    print("multi_processing ::: get_PAM_within_N_bp_of_PAM >>> ")
    result_list = []
    logic = Logic.Logics()

    for mut_arr in mut_list:
        chr_num = mut_arr[0]
        pos = int(mut_arr[1]) - 1
        ref_p_seq = mut_arr[3]
        alt_p_seq = mut_arr[4]

        ref_m_seq = ""
        alt_m_seq = ""
        try:
            ref_m_seq += logic.make_complement_string(ref_p_seq)
            if alt_p_seq == '.':
                alt_p_seq = ""
            else:
                alt_m_seq += logic.make_complement_string(alt_p_seq)
        except Exception as err:
            print("make_complement_string ::: ", err)
            print(ref_p_seq, " : ref_p_seq")
            print(alt_p_seq, " : alt_p_seq")
            print(str(mut_arr))

        seq_record = SeqIO.read(REF_DIR + "chr" + chr_num + ".fa", "fasta")
        p_seq = str(seq_record.seq).upper()
        m_seq = str(seq_record.seq.complement()).upper()

        logic.get_matched_pam_p_seq_list(result_list, p_seq, pos, PAM, WIN_PAM, SEQ_WIN_SIZE, "chr" + chr_num)
        logic.get_matched_pam_m_seq_list(result_list, m_seq, pos, PAM, WIN_PAM, SEQ_WIN_SIZE, "chr" + chr_num)

    print("DONE multi_processing ::: get_PAM_within_N_bp_of_PAM >>>")
    return result_list

def test():
    logic = Logic.Logics()
    util = Util.Utils()

    # CHROM	POS	ID	REF	ALT	mut_length	CLNVC	CLNSIG
    # POS - 1 = index of .fa sequence
    # [['1', '930188', '846933', 'G', 'A', '1', 'substitution', 'Uncertain_significance'],...]
    mut_list = util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t")[:20]

    result_list = get_PAM_within_N_bp_of_POS(mut_list)

    header = ['CHROM', 'PAM', str(SEQ_WIN_SIZE[0]) + ' + PAM + ' + str(SEQ_WIN_SIZE[1]), 'PAM_POS', 'STRAND']
    util.make_excel(WORK_DIR + "output/ClinVar_hg38_test_result", header, result_list)

if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    multi_processing()
    # test()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
