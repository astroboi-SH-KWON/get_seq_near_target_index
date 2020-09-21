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
    REF_DIR = "../hg38/"
else:
    # DEV
    REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"

PROJECT_NAME = WORK_DIR.split("/")[-2]
CDS_INFO = "hg38_refFlat_full.txt"
MUT_INFO = "20200909_ClinVar_hg38.txt"
STRT_CD_ARR = ['ATG']
END_CD_ARR = ['TGA', 'TAG', 'TAA']

PAM = 'NGG'
WIN_PAM = [10, 10]
SEQ_WIN_SIZE = [30, 30]


TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
############### end setting env #################

def multi_processing_1():
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

def test_get_PAM_within_N_bp_of_POS():
    logic = Logic.Logics()
    util = Util.Utils()

    # CHROM	POS	ID	REF	ALT	mut_length	CLNVC	CLNSIG
    # POS - 1 = index of .fa sequence
    # [['1', '930188', '846933', 'G', 'A', '1', 'substitution', 'Uncertain_significance'],...]
    mut_list = util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t")[:20]

    result_list = get_PAM_within_N_bp_of_POS(mut_list)

    header = ['CHROM', 'PAM', str(SEQ_WIN_SIZE[0]) + ' + PAM + ' + str(SEQ_WIN_SIZE[1]), 'PAM_POS', 'STRAND']
    util.make_excel(WORK_DIR + "output/ClinVar_hg38_test_result", header, result_list)

def make_filtered_hg38_refFlat():
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    # GeneSym   NMID    Chrom   Strand  Transcript_Start   End ORFStart    End #Exon   ExonS_list  ExonE_list
    # ['MIR6859-1', 'NR_106918', 'chr1', '-', '17368', '17436', '17436', '17436', '1', '17368,', '17436,']
    # ['WASH7P', 'NR_024540', 'chr1', '-', '14361', '29370', '29370', '29370', '11', '14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320,', '14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370,']
    cds_list = []
    if SYSTEM_NM == 'Linux':
        cds_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/" + CDS_INFO, "\t", 0))
    else:
        cds_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/" + CDS_INFO, "\t", 0)[:3000])

    NM_cds_list = logic_prep.filter_out_NON_NM_id_in_cds_list(cds_list)
    # filter_out_cds_wout_strt_cdn(NM_cds_list)

    splited_cds_list = np.array_split(NM_cds_list, MULTI_CNT)

    print("platform.system() : ", SYSTEM_NM)
    print("total cpu_count : ", str(TOTAL_CPU))
    print("will use : ", str(MULTI_CNT))
    pool = mp.Pool(processes=MULTI_CNT)

    pool_cds_idx_list = pool.map(filter_out_cds_wout_strt_cdn, splited_cds_list)
    result_list = logic_prep.merge_multi_list(pool_cds_idx_list)

    header = ['GeneSym', 'NMID', 'Chrom', 'Strand', 'Transcript_Start', 'End', 'ORFStart', 'End', '#Exon', 'ExonS_list',
              'ExonE_list', 'int(NMID)']
    util.make_excel(WORK_DIR + "output/filtered_hg38_refFlat", header, result_list)


def filter_out_cds_wout_strt_cdn(cds_list):
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()

    print("start filter_out_cds_wout_strt_cdn!!!!")
    result_list = []
    for cds_arr in cds_list:
        gene_sym = cds_arr[0]
        nm_id = cds_arr[1]
        chr_nm = cds_arr[2]
        strand = cds_arr[3]
        orf_strt_pos = int(cds_arr[6])
        orf_end_pos = int(cds_arr[7])

        p_seq, m_seq = util.read_file_by_biopython(REF_DIR + chr_nm + ".fa", "fasta")

        if strand == '+':
            strt_codon = p_seq[orf_strt_pos: orf_strt_pos + 3]
            if strt_codon in STRT_CD_ARR:
                end_codon = p_seq[orf_end_pos - 3: orf_end_pos]
                if end_codon in END_CD_ARR:
                    start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)

                    p_cds_seq = logic_prep.get_seq_by_idx_arr(p_seq, start_idx_arr, end_idx_arr)
                    if len(p_cds_seq) % 3 != 0:
                        continue
                    if logic.exist_another_orf_end_codon_in_cds_seq(p_cds_seq):
                        continue

                    tmp_arr = []
                    tmp_arr.extend(cds_arr)
                    result_list.append(tmp_arr)

        else:
            strt_codon = m_seq[orf_end_pos - 3: orf_end_pos][::-1]
            if strt_codon in STRT_CD_ARR:
                end_codon = m_seq[orf_strt_pos: orf_strt_pos + 3][::-1]
                if end_codon in END_CD_ARR:
                    start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)

                    m_cds_seq = logic_prep.get_seq_by_idx_arr(m_seq, start_idx_arr, end_idx_arr)
                    if len(m_cds_seq) % 3 != 0:
                        continue
                    if logic.exist_another_orf_end_codon_in_cds_seq(m_cds_seq, False):
                        continue

                    tmp_arr = []
                    tmp_arr.extend(cds_arr)
                    result_list.append(tmp_arr)

    print("DONE filter_out_cds_wout_strt_cdn!!!!")
    return result_list

def test():
    logic = Logic.Logics()
    test_list = [-141, 92, 182, 51, 125, 90, 186, 163, 116, 79, 500, 125, 111, 246]
    print(logic.is_all_bigger_than_cut_off(test_list))

if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # multi_processing_1()
    make_filtered_hg38_refFlat()
    # test()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
