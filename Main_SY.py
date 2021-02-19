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
# TYPE = 'mouse'
TYPE = 'human'
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

IN = 'input/'
OU = 'output/'
PROJECT_NAME = WORK_DIR.split("/")[-2]
MUT_INFO = "200907_Dominant filter.txt"
# FILTERED_CDS_INFO = "filtered_hg38_refFlat.txt"
# FILTERED_CDS_INFO = "filtered_CCDS.current.txt"  # 20201123 mouse
# FILTERED_CDS_INFO = "filtered_201130_CCDS_" + TYPE + "_current.txt"  # 20201130
FILTERED_CDS_INFO = "filtered_shortest_cdn_CCDS_" + TYPE + ".txt"  # 20201201
RESULT_FILE = "SY_Dominant_result_by_spacer.txt"

OTHOLOG = ['SaCas9', 'SaCas9_KKH', 'SaCas9_NNG', 'St1Cas9', 'Nm1Cas9', 'Nm2Cas9', 'CjCas9', 'SauriCas9', 'SauriCas9-KKH']
PAMS = ['NNGRRT', 'NNNRRT', 'NNG', 'NNRGAA', 'NNNNGATT', 'NNNNCC', 'NNNNRYAC', 'NNGG', 'NNRG']
SEQ_F_PAM = [43, 43, 43, 41, 45, 44, 44, 43, 43]
SEQ_B_PAM = [3, 3, 3, 3, 3, 3, 3, 3, 3]
WIN_SIZE = [60, 60]

RATIO = [0.05, 0.65]
ADJ_REF_IDX = -1
INIT = [OTHOLOG, PAMS, SEQ_F_PAM, SEQ_B_PAM, ADJ_REF_IDX]
# INIT_BY_PAM = [
#     ['SaCas9', 'NNGRRT', 43, 3, WIN_SIZE, RATIO]
#     , ['SaCas9_KKH', 'NNNRRT', 43, 3, WIN_SIZE, RATIO]
#     , ['SaCas9_NNG', 'NNG', 43, 3, WIN_SIZE, RATIO]
#     , ['St1Cas9', 'NNRGAA', 41, 3, WIN_SIZE, RATIO]
#     , ['Nm1Cas9', 'NNNNGATT', 45, 3, WIN_SIZE, RATIO]
#     , ['Nm2Cas9', 'NNNNCC', 44, 3, WIN_SIZE, RATIO]
#     , ['CjCas9', 'NNNNRYAC', 44, 3, WIN_SIZE, RATIO]
#     , ['SauriCas9', 'NNGG', 43, 3, WIN_SIZE, RATIO]
#     , ['SauriCas9-KKH', 'NNRG', 43, 3, WIN_SIZE, RATIO]
# ]

INIT_BY_PAM = [
['SaCas9', 'NNGRRT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9', 'NNGAAR', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9', 'NNGAAC', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9', 'NNGRGV', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9', 'NNGGAR', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-KKH', 'NNNAGT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-KKH', 'NNVGGT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-KKH', 'NNRAAT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-KKH', 'NNAGAT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-KKH', 'NNCRAT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-KKH', 'NNGGAT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-KKH', 'NNTGGT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-KKH', 'NNGARR', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-KKH', 'NNGAAC', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-KKH', 'NNTAAT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-KKH', 'NNGCGT', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-NNG', 'NNGH', 43, 3, WIN_SIZE, RATIO]
    , ['SaCas9-NNG', 'NNGG', 43, 3, WIN_SIZE, RATIO]
    , ['SauriCas9', 'NNGG', 43, 3, WIN_SIZE, RATIO]
    , ['SauriCas9', 'NNGA', 43, 3, WIN_SIZE, RATIO]
    , ['SauriCas9-KKH', 'NNRG', 43, 3, WIN_SIZE, RATIO]
    , ['SauriCas9-KKH', 'NNCG', 43, 3, WIN_SIZE, RATIO]
    , ['SauriCas9-KKH', 'NNVA', 43, 3, WIN_SIZE, RATIO]
    , ['St1Cas9', 'NNRGAA', 41, 3, WIN_SIZE, RATIO]
    , ['St1Cas9', 'NNAGGA', 41, 3, WIN_SIZE, RATIO]
    , ['St1Cas9', 'NNABCA', 41, 3, WIN_SIZE, RATIO]
    , ['St1Cas9', 'NNAWAA', 41, 3, WIN_SIZE, RATIO]
    , ['St1Cas9', 'NNGGGA', 41, 3, WIN_SIZE, RATIO]
    , ['St1Cas9', 'NNAGAB', 41, 3, WIN_SIZE, RATIO]
    , ['St1Cas9', 'NNCGAA', 41, 3, WIN_SIZE, RATIO]
    , ['Nm1Cas9', 'NNNNGATTD', 45, 3, WIN_SIZE, RATIO]
    , ['Nm1Cas9', 'NNNNGACTW', 45, 3, WIN_SIZE, RATIO]
    , ['Nm1Cas9', 'NNNNGYTTA', 45, 3, WIN_SIZE, RATIO]
    , ['Nm1Cas9', 'NNNNGYTTT', 45, 3, WIN_SIZE, RATIO]
    , ['Nm1Cas9', 'NNNNGATTC', 45, 3, WIN_SIZE, RATIO]
    , ['Nm1Cas9', 'NNNNGTCTA', 45, 3, WIN_SIZE, RATIO]
    , ['Nm1Cas9', 'NNNNGATAG', 45, 3, WIN_SIZE, RATIO]
    , ['Nm1Cas9', 'NNNNGAGTA', 45, 3, WIN_SIZE, RATIO]
    , ['Nm2Cas9',  'NNNNCCA', 44, 3, WIN_SIZE, RATIO]
    , ['Nm2Cas9',  'NNNNCCB', 44, 3, WIN_SIZE, RATIO]
    , ['CjCas9', 'NNNNACAC', 44, 3, WIN_SIZE, RATIO]
    , ['CjCas9', 'NNNNRTAC', 44, 3, WIN_SIZE, RATIO]
    , ['CjCas9', 'NNNNGCAC', 44, 3, WIN_SIZE, RATIO]
    , ['CjCas9', 'NNNNACAT', 44, 3, WIN_SIZE, RATIO]
    , ['CjCas9', 'NNNNGTAT', 44, 3, WIN_SIZE, RATIO]
]

# multi_processing_ClinVar_by_all_cds()
# ALL_CDS_INFO = "all_ccds_filtered_201130_CCDS_mouse_current.txt"  # mouse
ALL_CDS_INFO = "all_ccds_filtered_201130_CCDS_human_current.txt"  # human
ORI_CDS_INFO = "201130_CCDS_mouse_current.txt"  # mouse
# ORI_CDS_INFO = "201130_CCDS_human_current.txt"  # human
FILTERED_MUT_INFO = "ClinVar_dominant_mutation_on_CDS.txt"
INIT_FOR_CLINVAR = [
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

TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.8)
ADJ_REF_IDX = -1
############### end setting env #################


def get_cds_idx_list():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    # cds_dict = logic_prep.get_dict_from_list_by_ele_key(
    #     util.read_csv_ignore_N_line(WORK_DIR + "input/" + ALL_CDS_INFO, "\t"), 2)
    ccds_list = util.read_csv_ignore_N_line(WORK_DIR + "input/" + ORI_CDS_INFO, "\t", 1)
    ccds_list = logic_prep.filter_out_data_with_trgt_strng(ccds_list, '-', -2)
    ccds_hg38_form_list = logic_prep.transform_mouse_ccds_form_to_hg38_refFlat(ccds_list)
    # util.make_csv(WORK_DIR + 'output/ccds_hg38_form_list.txt', [], ccds_hg38_form_list, deli='\t')
    cds_dict = logic_prep.get_dict_from_list_by_ele_key(ccds_hg38_form_list, 2)

    cds_info = cds_dict['chr1']
    cds_idx_list = []
    for cds_arr in cds_info:
        start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)
        idx_list = logic_prep.get_idx_num_frm_strt_to_end_list(start_idx_arr, end_idx_arr)
        cds_idx_list.append(idx_list)
    print()
    print(logic.check_seq_in_cds(cds_idx_list, 97104905))
    util.make_csv(WORK_DIR + 'input/cds_idx_list_.txt', [], cds_idx_list)


def get_GRCh38_Regulatory_Build_regulatory_features_by_ClinVar_dominant_mutation_not_on_CDS():
    clin_var_not_cds_fl_nm = 'ClinVar_dominant_mutation_not_on_CDS.txt'
    GRCh38_features_fl_nm = 'homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff'

    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    clin_var_fl = util.read_csv_ignore_N_line(WORK_DIR + IN + clin_var_not_cds_fl_nm, '\t')
    GRCh38_features_fl = util.read_csv_ignore_N_line(WORK_DIR + IN + GRCh38_features_fl_nm, '\t', n_line=0)
    GRCh38_features_dict = logic_prep.get_dict_from_list_by_ele_key(GRCh38_features_fl, 0)

    result_dict = {}
    no_chr_key_list = []
    for clin_var_arr in clin_var_fl:
        tmp_clin_var_key = tuple(clin_var_arr[:5])
        chr_key = clin_var_arr[0]
        pos = int(clin_var_arr[1])
        if chr_key in GRCh38_features_dict:
            result_dict.update({tmp_clin_var_key: []})
            GRCh38_features_list = GRCh38_features_dict[chr_key]
            for GRCh38_features_arr in GRCh38_features_list:
                tmp_type = GRCh38_features_arr[2]
                tmp_info = GRCh38_features_arr[-1]
                st_idx = int(GRCh38_features_arr[3])
                en_idx = int(GRCh38_features_arr[4])
                if st_idx < pos < en_idx:
                    result_dict[tmp_clin_var_key].append([tmp_type, tmp_info])

        else:
            no_chr_key_list.append(tmp_clin_var_key)

    with open(WORK_DIR + OU + 'in_cds.txt', 'w') as in_cds_f:
        with open(WORK_DIR + OU + 'not_in_cds.txt', 'w') as not_cds_f:
            for f_key, val_list in result_dict.items():
                tmp_str = ""
                tmp_blnk = ""
                for tmp_f in f_key:
                    tmp_str = tmp_str + tmp_f + '\t'
                    tmp_blnk = tmp_blnk + '-' + '\t'

                if len(val_list) == 0:
                    not_cds_f.write(tmp_str[:-1] + '\n')
                else:
                    for idx in range(len(val_list)):
                        add_str = ""
                        for tm_f in val_list[idx]:
                            add_str = add_str + tm_f + '\t'

                        if idx == 0:
                            in_cds_f.write(tmp_str + add_str[:-1] + '\n')
                        else:
                            in_cds_f.write(tmp_blnk + add_str[:-1] + '\n')

    util.make_csv(WORK_DIR + OU + 'no_chr_key_list.txt', [], no_chr_key_list)


def make_filtered_out_ClinVar_pos_in_cds_or_not():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    cds_info = util.read_csv_ignore_N_line(WORK_DIR + "input/" + ALL_CDS_INFO, "\t")
    cds_dict_by_chr = {}
    for cds_arr in cds_info:
        chrom = cds_arr[2]
        start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)
        idx_list = logic_prep.get_idx_num_frm_strt_to_end_list(start_idx_arr, end_idx_arr)
        if chrom in cds_dict_by_chr:
            cds_dict_by_chr[chrom].append(idx_list)
        else:
            cds_dict_by_chr.update({chrom: [idx_list]})

    mut_dict = logic_prep.get_dict_from_list_by_ele_key(
        util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t"), 0)

    not_in_cds_list = []
    in_cds_list = []
    for chr_num, mut_list in mut_dict.items():
        cds_idx_list = cds_dict_by_chr['chr' + chr_num]
        for mut_arr in mut_list:
            pos = int(mut_arr[1]) + ADJ_REF_IDX
            tmp_id = int(mut_arr[2])
            if not logic.check_seq_in_cds(cds_idx_list, pos):
                not_in_cds_list.append(mut_arr)
            else:
                in_cds_list.append(mut_arr)
    print(len(not_in_cds_list))
    header = ['#CHROM',	'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    util.make_csv(WORK_DIR + '/input/ClinVar_dominant_mutation_on_CDS.txt', header, in_cds_list, deli='\t')
    util.make_csv(WORK_DIR + '/input/ClinVar_dominant_mutation_not_on_CDS.txt', header, not_in_cds_list, deli='\t')


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
        print("\nstart to make files for ", pam_nm, init[1] + "\n")
        util.make_csv(WORK_DIR + "output/cleavage_pos_in_shortest_cds_w_whole_" + TYPE + "_gene_" + pam_nm + '_' + init[1] + ".txt", header,
                      sorted_result_list, 0, "\t")
        # try:
        #     print(len(sorted_result_list))
        #     if pam_nm in ['SaCas9_NNG', 'SauriCas9-KKH']:
        #         continue
        #     util.make_excel(WORK_DIR + "output/cleavage_pos_in_shortest_cds_w_whole_" + TYPE + "_gene_" + pam_nm, header,
        #                     sorted_result_list)
        # except Exception as err:
        #     print('[ERROR] during making excel file for ' + pam_nm + '\n', str(err))


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
    # multi_processing_for_whole_pam_ClinVar()
    make_filtered_out_ClinVar_pos_in_cds_or_not()
    # get_cds_idx_list()
    # get_GRCh38_Regulatory_Build_regulatory_features_by_ClinVar_dominant_mutation_not_on_CDS()
    # multi_processing_ClinVar_by_all_cds()
    # split_multi_processing_for_whole_pam_by_PAM()
    # multi_processing_by_pam_w_whole_gene()
    # make_filtered_mouse_ccds_current_file()
    # make_filtered_ccds_current_file_by_shortest_cdn()  # 20201201
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
