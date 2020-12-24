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
TYPE = 'mouse'
# TYPE = 'human'
if SYSTEM_NM == 'Linux':
    # REAL
    # REF_DIR = "/media/backup/ref/hg38/"
    # REF_DIR = "/media/backup/ref/GRCm38_p6_mouse/"  # 20201123 mouse
    if TYPE == 'mouse':
        REF_DIR = "/media/backup/ref/Ensemble_GRCm38_p6/"  # 20201130 mouse
    elif TYPE == 'human':
        REF_DIR = "/media/backup/ref/Ensemble_GRCh38_p13/"  # 20201130 human
    else:
        import sys
        print("[ERROR] - no ref path")
        sys.exit()
    print('type : [', TYPE, ']')
else:
    # DEV
    REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"
    WORK_DIR = "D:/000_WORK/SeoSangYeon/20200907_ClinVar/WORK_DIR/"

PROJECT_NAME = WORK_DIR.split("/")[-2]
MUT_INFO = "200907_Dominant filter.txt"
# FILTERED_CDS_INFO = "filtered_hg38_refFlat.txt"
# FILTERED_CDS_INFO = "filtered_CCDS.current.txt"  # 20201123 mouse
# FILTERED_CDS_INFO = "filtered_201130_CCDS_" + TYPE + "_current.txt"  # 20201130
FILTERED_CDS_INFO = "filtered_shortest_cdn_CCDS_" + TYPE + ".txt"  # 20201201
RESULT_FILE = "SY_Dominant_result_by_spacer.txt"

# multi_processing_ClinVar_by_all_cds()
ALL_CDS_INFO = "all_ccds_filtered_201130_CCDS_" + TYPE + "_current.txt"  # 20201217

OTHOLOG = ['SaCas9', 'SaCas9_KKH', 'SaCas9_NNG', 'St1Cas9', 'Nm1Cas9', 'Nm2Cas9', 'CjCas9', 'SauriCas9', 'SauriCas9-KKH']
PAMS = ['NNGRRT', 'NNNRRT', 'NNG', 'NNRGAA', 'NNNNGATT', 'NNNNCC', 'NNNNRYAC', 'NNGG', 'NNRG']
SEQ_F_PAM = [43, 43, 43, 41, 45, 44, 44, 43, 43]
SEQ_B_PAM = [3, 3, 3, 3, 3, 3, 3, 3, 3]
WIN_SIZE = [60, 60]

RATIO = [0.05, 0.65]
ADJ_REF_IDX = -1
INIT = [OTHOLOG, PAMS, SEQ_F_PAM, SEQ_B_PAM, ADJ_REF_IDX]
INIT_BY_PAM = [
    ['SaCas9', 'NNGRRT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9_KKH', 'NNNRRT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9_NNG', 'NNG', 43, 3, WIN_SIZE, RATIO]
    , ['St1Cas9', 'NNRGAA', 41, 3, WIN_SIZE, RATIO]
    , ['Nm1Cas9', 'NNNNGATT', 45, 3, WIN_SIZE, RATIO]
    , ['Nm2Cas9', 'NNNNCC', 44, 3, WIN_SIZE, RATIO]
    , ['CjCas9', 'NNNNRYAC', 44, 3, WIN_SIZE, RATIO]
    , ['SauriCas9', 'NNGG', 43, 3, WIN_SIZE, RATIO]
    , ['SauriCas9-KKH', 'NNRG', 43, 3, WIN_SIZE, RATIO]
]

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)

############### end setting env #################

def multi_processing_ClinVar_by_all_cds():
    pass


def multi_processing_by_pam_w_whole_gene():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    file_nm_arr = ['chrX', 'chrY']
    if SYSTEM_NM == 'Linux':
        for f_num in range(1, 23):
            file_nm_arr.append("chr" + str(f_num))
    else:
        for f_num in range(20, 23):
            file_nm_arr.append("chr" + str(f_num))

    cds_info = util.read_csv_ignore_N_line(WORK_DIR + "input/" + FILTERED_CDS_INFO, "\t")

    cds_dict = logic_prep.get_dict_from_list_by_ele_key(cds_info, 2)

    for init in INIT_BY_PAM:
        pam_nm = init[0]
        len_f_pam = init[2]
        len_b_pam = init[3]
        init.extend([cds_dict])

        manager = mp.Manager()
        result_list = manager.list()
        jobs = []
        for key in file_nm_arr:
            if key not in cds_dict:
                continue
            proc = mp.Process(target=logic.get_seq_clvg_pos_in_cds_w_whole_gene, args=(REF_DIR, key, init, result_list))
            jobs.append(proc)
            proc.start()

        for ret_proc in jobs:
            ret_proc.join()

        header = ['chr_nm', 'gene_sym', 'nm_id', 'transcrpt_strand', 'pos_clvg', 'pam_strand', 'context (22 bp)',
                  'spacer (21 bp)', 'PAM', '3 bp', str(len_f_pam) + ' bp + PAM + ' + str(len_b_pam) + ' bp',
                  'clvg_site_ratio']

        sorted_result_list = logic_prep.sort_list_by_ele(result_list, 0, False)
        # clear ListProxy
        result_list[:] = []
        print("\nstart to make files for ", pam_nm, "\n")
        util.make_csv(WORK_DIR + "output/cleavage_pos_in_shortest_cds_w_whole_" + TYPE + "_gene_" + pam_nm + ".txt", header,
                      sorted_result_list, 0, "\t")
        try:
            print(len(sorted_result_list))
            if pam_nm in ['SaCas9_NNG', 'SauriCas9-KKH']:
                continue
            util.make_excel(WORK_DIR + "output/cleavage_pos_in_shortest_cds_w_whole_" + TYPE + "_gene_" + pam_nm, header,
                            sorted_result_list)
        except Exception as err:
            print('[ERROR] during making excel file for ' + pam_nm + '\n', str(err))


def multi_processing_by_pam_w_mut():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    mut_list = []
    file_nm_arr = ['chrX', 'chrY']
    if SYSTEM_NM == 'Linux':
        mut_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t"))
        for f_num in range(1, 23):
            file_nm_arr.append("chr" + str(f_num))
    else:
        # 20, 21, 22, X
        mut_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t")[20000:])
        for f_num in range(20, 23):
            file_nm_arr.append("chr" + str(f_num))

    cds_info = util.read_csv_ignore_N_line(WORK_DIR + "input/" + FILTERED_CDS_INFO, "\t")

    mut_dict = logic_prep.get_dict_from_list_by_ele_key(mut_list, 0)
    cds_dict = logic_prep.get_dict_from_list_by_ele_key(cds_info, 2)

    for init in INIT_BY_PAM:
        pam_nm = init[0]
        len_f_pam = init[2]
        len_b_pam = init[3]
        init.extend([cds_dict, mut_dict])

        manager = mp.Manager()
        result_list = manager.list()
        jobs = []
        for key in file_nm_arr:
            if key not in cds_dict:
                continue
            if key.replace("chr", "") not in mut_dict:
                continue
            proc = mp.Process(target=logic.get_seq_clvg_pos_in_cds, args=(REF_DIR, key, init, result_list))
            jobs.append(proc)
            proc.start()

        for ret_proc in jobs:
            ret_proc.join()

        header = ['chr_nm', 'gene_sym', 'nm_id', 'pos_clvg', 'strand', '22 bp', 'spacer', 'PAM', '3bp',
                  str(len_f_pam) + 'bp + PAM + ' + str(len_b_pam) + 'bp', 'clvg_site_ratio']
        sorted_result_list = logic_prep.sort_list_by_ele(result_list, 0, False)
        print("\nstart to make files for ", pam_nm)
        util.make_csv(WORK_DIR + "output/cleavage_pos_in_cds_w_mut_" + pam_nm + ".txt", header, sorted_result_list, 0,
                      "\t")
        util.make_excel(WORK_DIR + "output/cleavage_pos_in_cds_w_mut_" + pam_nm, header, sorted_result_list)


def multi_processing_for_whole_pam_ClinVar():
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    # POS - 1 = index of .fa sequence
    # [['1', '1338000', '208047', 'CT', 'C', '.', '.', '"ALLELEID=204306;CLNDISDB=MONDO:MONDO:0014591,MedGen:C4225363,OMIM:616331;CLNDN=Robinow_syndrome,_autosomal_dominant_2;CLNHGVS=NC_000001.11:g.1338001del;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Pathogenic;CLNVC=Deletion;CLNVCSO=SO:0000159;GENEINFO=DVL1:1855;MC=SO:0001589|frameshift_variant;ORIGIN=33;RS=797044837"'],...]
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

    pool_list = pool.map(get_seq_by_pam_after_mut, splited_mut_list)

    result_list = logic_prep.merge_multi_list(pool_list)

    header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
              'P_REF_SEQ_[' + str(WIN_SIZE[0]) + '], M_REF_SEQ_[' + str(WIN_SIZE[1]) + ']']
    for pam_nm in OTHOLOG:
        for strand in ['+', '-']:
            header.append(pam_nm + strand)
    util.make_excel(WORK_DIR + "output/SY_Dominant_result_by_spacer", header, result_list)


def get_seq_by_pam_after_mut(mut_list):
    print("multi_processing ::: get_seq_by_pam_after_mut >>> ")
    result_list = []
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()

    win_arr = WIN_SIZE
    init = INIT

    pam_arr = init[1]
    len_f_pam_arr = init[2]
    len_b_pam_arr = init[3]

    for mut_arr in mut_list:
        tmp_arr = []
        tmp_arr.extend(mut_arr)
        chr_num = mut_arr[0]
        pos = int(mut_arr[1]) + ADJ_REF_IDX
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

        ori_win_flag = True

        for idx in range(len(pam_arr)):
            pam = pam_arr[idx]
            len_f_pam = len_f_pam_arr[idx]
            len_b_pam = len_b_pam_arr[idx]

            ref_p_dict, p_ori_win_seq = logic.get_matched_pam_p_seq_dict(p_seq, pos, win_arr, ref_p_seq, pam,
                                                                         len_f_pam, len_b_pam)
            ref_m_dict, m_ori_win_seq = logic.get_matched_pam_m_seq_dict(m_seq, pos, win_arr, ref_m_seq, pam,
                                                                         len_f_pam, len_b_pam)

            mut_p_dict, _ = logic.get_matched_pam_p_seq_dict(p_seq, pos, win_arr, alt_p_seq, pam, len_f_pam,
                                                             len_b_pam)
            mut_m_dict, _ = logic.get_matched_pam_m_seq_dict(m_seq, pos, win_arr, alt_m_seq, pam, len_f_pam,
                                                             len_b_pam)

            logic.remove_dict_val_by_key(mut_p_dict, ref_p_dict.keys())
            logic.remove_dict_val_by_key(mut_m_dict, ref_m_dict.keys())

            if ori_win_flag:
                tmp_arr.append(p_ori_win_seq + " , " + m_ori_win_seq)
                ori_win_flag = False

            logic_prep.add_result_seq_to_arr(tmp_arr, mut_p_dict)
            logic_prep.add_result_seq_to_arr(tmp_arr, mut_m_dict)

        result_list.append(tmp_arr)
    print("DONE multi_processing ::: get_seq_by_pam_after_mut >>>")
    return result_list


def split_multi_processing_for_whole_pam_by_PAM():
    util = Util.Utils()
    result_list = util.read_csv_ignore_N_line(WORK_DIR + "input/" + RESULT_FILE, "\t")
    result_dict = {
        'SaCas9': []
        , 'SaCas9_KKH': []
        , 'SaCas9_NNG': []
        , 'St1Cas9': []
        , 'Nm1Cas9': []
        , 'Nm2Cas9': []
        , 'CjCas9': []
    }
    key_arr = ['SaCas9', 'SaCas9',
              'SaCas9_KKH', 'SaCas9_KKH', 'SaCas9_NNG', 'SaCas9_NNG', 'St1Cas9', 'St1Cas9', 'Nm1Cas9',
              'Nm1Cas9', 'Nm2Cas9', 'Nm2Cas9', 'CjCas9', 'CjCas9']
    for val_arr in result_list:
        ori_p = val_arr[8].replace('"', '').split(",")[0]
        ori_m = val_arr[8].replace('"', '').split(",")[1]

        for i in range(9, 23):
            tmp_str = val_arr[i]
            if tmp_str != '':
                for val in tmp_str.replace('"', '').split(",")[:-1]:
                    tmp_arr = val_arr[:8]
                    if i % 2 == 0:
                        tmp_arr.append("-")
                        tmp_arr.append(ori_m)
                    else:
                        tmp_arr.append("+")
                        tmp_arr.append(ori_p)
                    tmp_arr.append(val)
                    result_dict[key_arr[i - 9]].append(tmp_arr)

    for key, val in result_dict.items():
        header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'strand',
                  '[' + str(WIN_SIZE[0]) + ']_REF_SEQ_[' + str(WIN_SIZE[1]) + ']', 'seq']

        util.make_csv(WORK_DIR + "output/SY_Dominant_result_by_spacer_" + key + ".txt", header, val, 0, "\t")
        util.make_excel(WORK_DIR + "output/SY_Dominant_result_by_spacer_" + key, header, val)

def main():
    logic = Logic.Logics()
    util = Util.Utils()

    # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    # POS - 1 = index of .fa sequence
    # [['1', '1338000', '208047', 'CT', 'C', '.', '.', '"ALLELEID=204306;CLNDISDB=MONDO:MONDO:0014591,MedGen:C4225363,OMIM:616331;CLNDN=Robinow_syndrome,_autosomal_dominant_2;CLNHGVS=NC_000001.11:g.1338001del;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Pathogenic;CLNVC=Deletion;CLNVCSO=SO:0000159;GENEINFO=DVL1:1855;MC=SO:0001589|frameshift_variant;ORIGIN=33;RS=797044837"'],...]
    mut_list = util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t")

    logic.get_seq_by_pam_after_mut(REF_DIR, mut_list, WIN_SIZE, INIT)

    header = ['#CHROM',	'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'P_REF_SEQ_[' + str(WIN_SIZE[0]) + '], M_REF_SEQ_[' + str(WIN_SIZE[1]), 'SaCas9+', 'SaCas9-', 'SaCas9_KKH+', 'SaCas9_KKH-', 'SaCas9_NNG+', 'SaCas9_NNG-', 'St1Cas9+', 'St1Cas9-', 'Nm1Cas9+', 'Nm1Cas9-', 'Nm2Cas9+', 'Nm2Cas9-', 'CjCas9+', 'CjCas9-']
    util.make_excel(WORK_DIR + "output/SY_Dominant_test_result_by_spacer", header, mut_list)


def make_filtered_mouse_ccds_current_file():
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    ccds_list = []
    if SYSTEM_NM == 'Linux':
        ccds_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/201130_CCDS_human_current.txt", "\t", 0))
    else:
        # ccds_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/CCDS.current.txt", "\t", 0)[:3000])
        ccds_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/201130_CCDS_human_current.txt", "\t", 0))

    # st plan A : filter out non Public, non Identical
    ccds_list = logic_prep.get_data_with_trgt_strng(ccds_list, 'Public', 5)
    ccds_list = logic_prep.get_data_with_trgt_strng(ccds_list, 'Identical', -1)
    # get the highest num of ccds_id in each gene
    ccds_list = logic_prep.get_highest_ccds_id_among_same_gen_id(ccds_list)
    # en plan A

    ccds_hg38_form_list = logic_prep.transform_mouse_ccds_form_to_hg38_refFlat(ccds_list)

    header = ['GeneSym', 'NMID', 'Chrom', 'Strand', 'Transcript_Start', 'End', 'ORFStart', 'End', '#Exon', 'ExonS_list',
              'ExonE_list']

    try:
        os.remove(WORK_DIR + "input/filtered_CCDS.current.txt")
    except Exception as err:
        print('os.remove(WORK_DIR + "input/filtered_CCDS.current.txt") : ', str(err))
    util.make_csv(WORK_DIR + "input/filtered_201130_CCDS_human_current.txt", header, ccds_hg38_form_list, 0, "\t")


def make_filtered_ccds_current_file_by_shortest_cdn():
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    ccds_list = []
    if SYSTEM_NM == 'Linux':
        ccds_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/201130_CCDS_" + TYPE + "_current.txt", "\t", 0))
    else:
        # ccds_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/CCDS.current.txt", "\t", 0)[:3000])
        ccds_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/201130_CCDS_" + TYPE + "_current.txt", "\t", 0))

    # st plan A : filter out non Public, non Identical
    ccds_list = logic_prep.get_data_with_trgt_strng(ccds_list, 'Public', 5)
    ccds_list = logic_prep.get_data_with_trgt_strng(ccds_list, 'Identical', -1)

    ccds_hg38_form_list = logic_prep.transform_mouse_ccds_form_to_hg38_refFlat(ccds_list)

    filted_ccds_list = logic_prep.get_shortest_cdn_among_same_gen_id(ccds_hg38_form_list)  # 20201201
    # en plan A

    header = ['GeneSym', 'NMID', 'Chrom', 'Strand', 'Transcript_Start', 'End', 'ORFStart', 'End', '#Exon', 'ExonS_list',
              'ExonE_list']

    try:
        os.remove(WORK_DIR + "input/filtered_shortest_cdn_CCDS_" + TYPE + ".txt")
    except Exception as err:
        print('os.remove(WORK_DIR + "input/filtered_CCDS.current.txt") : ', str(err))
    util.make_csv(WORK_DIR + "input/filtered_shortest_cdn_CCDS_" + TYPE + ".txt", header, filted_ccds_list, 0, "\t")


def check_seq_idx():
    util = Util.Utils()

    # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    # POS - 1 = index of .fa sequence
    # [['1', '1338000', '208047', 'CT', 'C', '.', '.', '"ALLELEID=204306;CLNDISDB=MONDO:MONDO:0014591,MedGen:C4225363,OMIM:616331;CLNDN=Robinow_syndrome,_autosomal_dominant_2;CLNHGVS=NC_000001.11:g.1338001del;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Pathogenic;CLNVC=Deletion;CLNVCSO=SO:0000159;GENEINFO=DVL1:1855;MC=SO:0001589|frameshift_variant;ORIGIN=33;RS=797044837"'],...]
    mut_list = []
    if SYSTEM_NM == 'Linux':
        mut_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t"))
    else:
        mut_list.extend(util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t")[:300])

    for mut_arr in mut_list:
        chr_num = mut_arr[0]
        pos = int(mut_arr[1]) + ADJ_REF_IDX
        ref_seq = mut_arr[3]
        seq_record = SeqIO.read(REF_DIR + "chr" + chr_num + ".fa", "fasta")
        print(ref_seq)
        print(str(seq_record.seq)[pos: pos + len(ref_seq)])
        print()



if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # check_seq_idx()
    multi_processing_for_whole_pam_ClinVar()
    multi_processing_ClinVar_by_all_cds()
    # split_multi_processing_for_whole_pam_by_PAM()
    # multi_processing_by_pam_w_whole_gene()
    # make_filtered_mouse_ccds_current_file()
    # make_filtered_ccds_current_file_by_shortest_cdn()  # 20201201
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
