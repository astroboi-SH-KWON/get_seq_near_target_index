import time
import os
from Bio import SeqIO
import multiprocessing as mp
import numpy as np
import platform

import Util
import Logic
import LogicPrep
############### st env ################
WORK_DIR = os.getcwd() + "/"
SYSTEM_NM = platform.system()
if SYSTEM_NM == 'Linux':
    # REAL
    # REF_DIR = "/media/backup/ref/hg38/"
    # REF_DIR = "/media/backup/ref/GRCm38_p6_mouse/"  # 20201123 mouse
    pass
else:
    # DEV
    WORK_DIR = "D:/000_WORK/ChangYoojin/20201207/WORK_DIR/"


IN = 'input/'
OU = 'output/'

SEQ_FL = 'Agt_rattus_norvegicus.txt'

PROJECT_NAME = WORK_DIR.split("/")[-2]

INIT = [
    ['SaCas9', 22, 21, 'NNGRRT', 3]
    , ['SaCas9-KKH', 22, 21, 'NNNRRT', 3]
    , ['SaCas9-NNG', 22, 21, 'NNG', 3]
    , ['SauriCas9', 22, 21, 'NNGG', 3]
    , ['SauriCas9-KKH', 22, 21, 'NNRG', 3]
    , ['St1Cas9', 22, 19, 'NNRGAA', 3]
    , ['Nm1Cas9', 22, 23, 'NNNNGATT', 3]
    , ['Nm2Cas9', 22, 22, 'NNNNCC', 3]
    , ['CjCas9', 22, 22, 'NNNNRYAC', 3]
]
############### en env ################


def main():
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()
    util = Util.Utils()

    result_list = []
    with open(WORK_DIR + IN + SEQ_FL) as f:
        p_seq_str = f.readline().upper()
        m_seq_str = logic.make_complement_string(p_seq_str)

        for idx in range(len(p_seq_str)):
            for pam_inf_arr in INIT:
                p_pam = pam_inf_arr[3]
                m_pam = p_pam[::-1]
                len_pam = len(p_pam)

                pam_fr_p_seq = p_seq_str[idx: idx + len_pam]
                pam_fr_m_seq = m_seq_str[idx: idx + len_pam]

                if len(pam_fr_p_seq) < len_pam or len(pam_fr_m_seq) < len_pam:
                    continue

                if logic.match(0, pam_fr_p_seq, p_pam):
                    pam_nm = pam_inf_arr[0]
                    len_f_pam = pam_inf_arr[1]
                    len_guide = pam_inf_arr[2]
                    len_b_pam = pam_inf_arr[4]
                    p_guide_seq = p_seq_str[idx - len_guide: idx]
                    if len(p_guide_seq) == len_guide:
                        result_list.append([pam_nm, p_guide_seq, pam_fr_p_seq, '+'])

                if logic.match(0, pam_fr_m_seq, m_pam):
                    pam_nm = pam_inf_arr[0]
                    len_f_pam = pam_inf_arr[1]
                    len_guide = pam_inf_arr[2]
                    len_b_pam = pam_inf_arr[4]
                    m_guide_seq = m_seq_str[idx + len_pam: idx + len_pam + len_guide][::-1]
                    if len(m_guide_seq) == len_guide:
                        result_list.append([pam_nm, m_guide_seq, pam_fr_m_seq[::-1], '-'])

    sorted_result_list = logic_prep.sort_list_by_ele(result_list, 0)
    util.make_csv(WORK_DIR + OU + 'result.txt', ['pam_name', 'guide_seq', 'pam_seq', 'strand'], sorted_result_list, deli='\t')




if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    main()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))