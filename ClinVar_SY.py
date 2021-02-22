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
    REF_DIR = "/media/backup/ref/Ensemble_GRCh38_p13/"  # 20201130 human
else:
    # DEV
    REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"
    WORK_DIR = "D:/000_WORK/SeoSangYeon/20200907_ClinVar/WORK_DIR/"

IN = 'input/'
OU = 'output/'
PROJECT_NAME = WORK_DIR.split("/")[-2]

ADJ_REF_IDX = -1
ALL_CDS_INFO = "all_ccds_filtered_201130_CCDS_human_current.txt"  # human
ORI_CDS_INFO = "201130_CCDS_human_current.txt"  # human
MUT_ON_CDS_INFO = "ClinVar_dominant_mutation_on_CDS.txt"

INIT_FOR_CLINVAR = [
    ['SaCas9', 22, 21, ['NNGRRT', 'NNGAAR', 'NNGAAC', 'NNGRGV', 'NNGGAR'], 3]
    , ['SaCas9-KKH', 22, 21, ['NNNAGT', 'NNVGGT', 'NNRAAT', 'NNAGAT', 'NNCRAT', 'NNGGAT', 'NNTGGT', 'NNGARR', 'NNGAAC', 'NNTAAT', 'NNGCGT'], 3]
    , ['SaCas9-NNG', 22, 21, ['NNGH', 'NNGG'], 3]
    , ['SauriCas9', 22, 21, ['NNGG', 'NNGA'], 3]
    , ['SauriCas9-KKH', 22, 21, ['NNRG', 'NNCG', 'NNVA'], 3]
    , ['St1Cas9', 22, 19, ['NNRGAA', 'NNAGGA', 'NNABCA', 'NNAWAA', 'NNGGGA', 'NNAGAB', 'NNCGAA'], 3]
    , ['Nm1Cas9', 22, 23, ['NNNNGATTD', 'NNNNGACTW', 'NNNNGYTTA', 'NNNNGYTTT', 'NNNNGATTC', 'NNNNGTCTA', 'NNNNGATAG', 'NNNNGAGTA'], 3]
    , ['Nm2Cas9', 22, 22, ['NNNNCCA', 'NNNNCCB'], 3]
    , ['CjCas9', 22, 22, ['NNNNACAC', 'NNNNRTAC', 'NNNNGCAC', 'NNNNACAT', 'NNNNGTAT'], 3]
]

WIN_SIZE = [60, 60]
TOTAL_CPU = mp.cpu_count()
MULTI_CNT = int(TOTAL_CPU*0.9)

all_cds_info = Util.Utils().read_csv_ignore_N_line(WORK_DIR + IN + ALL_CDS_INFO, deli_str='\t')
CDS_DICT = LogicPrep.LogicPreps().get_dict_from_list_by_ele_key(all_cds_info, 2)
all_cds_info.clear()
############### end setting env #################


def get_guide_set_from_ref(p_sq, m_sq_no_rvrsed, gene_list, cds_dict, init_arr):
    logic = Logic.Logics()

    result_set = set()

    # len_f_guide = init_arr[1]
    len_guide = init_arr[2]
    pam_rule_arr = init_arr[3]
    # len_b_pam = init_arr[4]

    for gen_nm in gene_list:
        cds_idx_list_arr = cds_dict[gen_nm]
        for cds_idx_list in cds_idx_list_arr:
            for clv_idx in cds_idx_list:

                for pam_rule in pam_rule_arr:
                    len_pam = len(pam_rule)

                    pam_fr_p_sq = p_sq[clv_idx + 3: clv_idx + 3 + len_pam]
                    pam_fr_m_sq = m_sq_no_rvrsed[clv_idx - 3 - len_pam: clv_idx - 3]

                    if logic.match_SY(0, pam_fr_p_sq, pam_rule):
                        guide_seq = p_sq[clv_idx + 3 - len_guide: clv_idx + 3]
                        result_set.add(guide_seq)

                    if logic.match_SY(0, pam_fr_m_sq, pam_rule[::-1]):
                        guide_seq = m_sq_no_rvrsed[clv_idx - 3: clv_idx - 3 + len_guide]
                        result_set.add(guide_seq[::-1])

    return result_set


def get_guide_dict_from_alt(p_sq, init_arr, mut_arr, gene_list, cds_dict):
    logic = Logic.Logics()
    result_dict = {}

    len_f_guide = init_arr[1]
    len_guide = init_arr[2]
    pam_rule_arr = init_arr[3]
    len_b_pam = init_arr[4]

    pos = int(mut_arr[1]) + ADJ_REF_IDX
    ref_p_seq = mut_arr[3]
    len_ref = len(ref_p_seq)
    alt_p_seq = mut_arr[4]
    # alt_p_seq == '.' ==> ''
    if alt_p_seq == '.':
        alt_p_seq = ''
    len_alt = len(alt_p_seq)

    # 2.3. ClinVar의 정보 (#CHROM, POS)를 기준으로 genome 상에서 mutation position의 주변 sequence를 가져옴 (60 bp + ref + 60 bp)
    p_f_win = p_sq[pos - WIN_SIZE[0]: pos]
    p_b_win = p_sq[pos + len_ref: pos + len_ref + WIN_SIZE[1]]
    p_win_seq_w_ref = p_f_win + ref_p_seq + p_b_win
    p_win_seq_w_alt = p_f_win + alt_p_seq + p_b_win

    for pam_rule in pam_rule_arr:
        len_pam = len(pam_rule)
        first_pam_pos = len(p_f_win) - len_pam + 1
        p_clvg_first_pos = pos - len_pam + 1 - 3
        m_clvg_first_pos = pos + 1 + 3

        # check PAM in only pos_seq part
        for i in range(len_pam + len_alt - 1):
            p_pam_seq = p_win_seq_w_alt[i + first_pam_pos: i + first_pam_pos + len_pam]
            m_pam_seq = logic.make_complement_string(p_pam_seq)

            p_clvg_pos = p_clvg_first_pos + i
            m_clvg_pos = m_clvg_first_pos + i

            if logic.match_SY(0, p_pam_seq, pam_rule) and logic.is_on_cds(p_clvg_pos, gene_list, cds_dict):
                f_guide = p_win_seq_w_alt[i + first_pam_pos - len_guide - len_f_guide: i + first_pam_pos - len_guide]
                guide_seq = p_win_seq_w_alt[i + first_pam_pos - len_guide: i + first_pam_pos]
                b_pam = p_win_seq_w_alt[i + first_pam_pos + len_pam: i + first_pam_pos + len_pam + len_b_pam]

                # STRAND, REF (2.3에서의 sequence), Guide context (22 nt + guide RNA + PAM + 3nt)
                tmp_arr = ['+', p_win_seq_w_ref, f_guide + guide_seq + p_pam_seq + b_pam]
                if guide_seq in result_dict:
                    result_dict[guide_seq].append(tmp_arr)
                else:
                    result_dict.update({guide_seq: [tmp_arr]})

            if logic.match_SY(0, m_pam_seq, pam_rule[::-1]) and logic.is_on_cds(m_clvg_pos, gene_list, cds_dict):
                m_win_seq_w_alt = logic.make_complement_string(p_win_seq_w_alt)
                f_guide = m_win_seq_w_alt[
                          i + first_pam_pos + len_pam + len_guide: i + first_pam_pos + len_pam + len_guide + len_f_guide]
                guide_seq = m_win_seq_w_alt[i + first_pam_pos + len_pam: i + first_pam_pos + len_pam + len_guide]
                b_pam = m_win_seq_w_alt[i + first_pam_pos - len_b_pam: i + first_pam_pos]

                # STRAND, REF (2.3에서의 sequence), Guide context (22 nt + guide RNA + PAM + 3nt)
                tmp_arr = ['-', logic.make_complement_string(p_win_seq_w_ref)[::-1],
                           (b_pam + m_pam_seq + guide_seq + f_guide)[::-1]]
                if guide_seq[::-1] in result_dict:
                    result_dict[guide_seq[::-1]].append(tmp_arr)
                else:
                    result_dict.update({guide_seq[::-1]: [tmp_arr]})

    return result_dict


def get_seq_by_pam_after_mut_ClinVar(mut_list):
    print("multi_processing ::: get_seq_by_pam_after_mut_ClinVar >>> ")
    result_list = []
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()

    # make mut_list to mut_dict by #CHROM (#CHROM doesn't have 'chr')
    mut_dict = logic_prep.get_dict_from_list_by_ele_key(mut_list, 0)

    for chr_key, tmp_mut_list in mut_dict.items():
        seq_record = SeqIO.read(REF_DIR + "chr" + chr_key + ".fa", "fasta")
        p_seq = str(seq_record.seq).upper()
        m_seq = str(seq_record.seq.complement()).upper()

        # GeneSym	NMID	Chrom	Strand	Transcript_Start	End	ORFStart	End	#Exon	ExonS_list	ExonE_list
        # Xkr4	CCDS14803.1	chr1	-	3216021	3671347	3216021	3671347	3	3216021,3421701,3670551,	3216967,3421900,3671347,
        cds_info = CDS_DICT["chr" + chr_key]
        cds_dict_by_GeneSym = logic.get_cds_idx_arr_dict_by_GeneSym(cds_info)

        for init_arr in INIT_FOR_CLINVAR:
            pam_nm = init_arr[0]

            # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            # POS - 1 = index of .fa sequence
            # [['1', '1338000', '208047', 'CT', 'C',
            for mut_arr in tmp_mut_list:
                pos = int(mut_arr[1]) + ADJ_REF_IDX
                ref_p_seq = mut_arr[3]
                # alt_p_seq = mut_arr[4]

                # 2.3.2. 이 때, ClinVar에서 보여주는 position 정보에 정확히 ref sequence가 있는지 검증하는 과정 필요
                if p_seq[pos: pos + len(ref_p_seq)] != ref_p_seq:
                    print('2.3.2. 이 때, ClinVar에서 보여주는 position 정보에 정확히 ref sequence가 있는지 검증하는 과정 필요\n', ref_p_seq,
                          ': ref_seq from ClinVar\n', p_seq[pos: pos + len(ref_p_seq)], ': ref_seq from .fa\n',
                          'mut_arr :', str(mut_arr), '\n\n')
                    continue

                # 2.1.2. 해당 filtered ClinVar 정보를 기준으로 각각의 mutation (ID)이 포함되는 CCDS id를 기준으로 해당 gene을 찾음
                GeneSym_list = []
                for cds_arr in cds_info:
                    GeneSym = cds_arr[0]
                    if GeneSym in GeneSym_list:
                        continue
                    start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)
                    for i in range(len(start_idx_arr)):
                        if start_idx_arr[i] <= pos <= end_idx_arr[i]:
                            GeneSym_list.append(GeneSym)
                            break
                if len(GeneSym_list) == 0:break

                ref_guide_set = get_guide_set_from_ref(p_seq, m_seq, GeneSym_list, cds_dict_by_GeneSym, init_arr)
                alt_guide_dict = get_guide_dict_from_alt(p_seq, init_arr, mut_arr, GeneSym_list, cds_dict_by_GeneSym)

                # 2.5. 2.4에서 가져온 sequence 중 guide를 기준 (주변의 sequence나 PAM은 포함하면 안됨)으로 2.1.에서 가져온 guide와 중복되는 것들은 제거하여 최종의 output을 만듦
                for guide_key, val_list in alt_guide_dict.items():
                    if guide_key not in ref_guide_set:
                        for val_arr in val_list:
                            # filter out if alt_guide in ref_seq
                            if logic.is_alt_guide_in_ref_seq(init_arr, val_arr):
                                break

                            tm_arr = []
                            tm_arr.extend(mut_arr)
                            tm_arr.extend(val_arr)
                            tm_arr.append(pam_nm)
                            result_list.append(tm_arr)
    print("DONE ::: get_seq_by_pam_after_mut_ClinVar >>> ")
    return result_list


def multi_processing_ClinVar_by_all_cds():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    clinvar_mut_info = util.read_csv_ignore_N_line(WORK_DIR + IN + MUT_ON_CDS_INFO, deli_str='\t')

    splited_clinvar_mut_info = np.array_split(clinvar_mut_info, MULTI_CNT)

    print("platform.system() : ", SYSTEM_NM)
    print("total cpu_count : ", str(TOTAL_CPU))
    print("will use : ", str(MULTI_CNT))
    pool = mp.Pool(processes=MULTI_CNT)

    pool_list = pool.map(get_seq_by_pam_after_mut_ClinVar, splited_clinvar_mut_info)

    result_list = logic_prep.merge_multi_list(pool_list)
    pool.close()
    pool_list[:] = []

    result_dict = logic_prep.get_dict_from_list_by_ele_key(result_list, -1)
    result_list.clear()

    header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'STRAND', 'REF (2.3에서의 sequence)', 'Guide context (22 nt + guide RNA + PAM + 3nt)']
    for pam_nm_key, rslt_list in result_dict.items():
        # remove pam_nm column
        fltd_rslt_list = [tmp_arr[:-1] for tmp_arr in rslt_list]
        print('make result file for', pam_nm_key)
        util.make_csv(WORK_DIR + OU + pam_nm_key + '_result.txt', header, fltd_rslt_list, deli='\t')


def get_guide_dict_from_ref(p_sq, m_sq_no_rvrsed, gene_list, cds_dict, init_arr):
    logic = Logic.Logics()

    result_dict = {}

    len_f_guide = init_arr[1]
    len_guide = init_arr[2]
    pam_rule = init_arr[3]
    len_pam = len(pam_rule)
    len_b_pam = init_arr[4]

    for gen_nm in gene_list:
        cds_idx_list_arr = cds_dict[gen_nm]
        for cds_idx_list in cds_idx_list_arr:
            for clv_idx in cds_idx_list:
                pam_fr_p_sq = p_sq[clv_idx + 3: clv_idx + 3 + len_pam]
                pam_fr_m_sq = m_sq_no_rvrsed[clv_idx - 3 - len_pam: clv_idx - 3]

                if logic.match_SY(0, pam_fr_p_sq, pam_rule):
                    f_guide = p_sq[clv_idx + 3 - len_guide - len_f_guide: clv_idx + 3 - len_guide]
                    guide_seq = p_sq[clv_idx + 3 - len_guide: clv_idx + 3]
                    b_pam = p_sq[clv_idx + 3 + len_pam: clv_idx + 3 + len_pam + len_b_pam]

                    # STRAND, Guide context (22 nt + guide RNA + PAM + 3nt)
                    tmp_arr = ['+', f_guide + guide_seq + pam_fr_p_sq + b_pam]
                    if guide_seq in result_dict:
                        result_dict[guide_seq].append(tmp_arr)
                    else:
                        result_dict.update({guide_seq: [tmp_arr]})

                if logic.match_SY(0, pam_fr_m_sq, pam_rule[::-1]):
                    f_guide = m_sq_no_rvrsed[clv_idx - 3 + len_guide: clv_idx - 3 + len_guide + len_f_guide]
                    guide_seq = m_sq_no_rvrsed[clv_idx - 3: clv_idx - 3 + len_guide]
                    b_pam = m_sq_no_rvrsed[clv_idx - 3 - len_pam - len_b_pam: clv_idx - 3 - len_pam]

                    # STRAND, Guide context (22 nt + guide RNA + PAM + 3nt)
                    tmp_arr = ['-', (b_pam + pam_fr_p_sq + guide_seq + f_guide)[::-1]]
                    if guide_seq in result_dict:
                        result_dict[guide_seq].append(tmp_arr)
                    else:
                        result_dict.update({guide_seq: [tmp_arr]})

    return result_dict


def test():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()

    # pos_arr = [26766456, 26772501]
    pos_arr = [33441359]
    chr_key = '6'

    cds_info = CDS_DICT["chr" + chr_key]
    cds_dict_by_GeneSym = logic.get_cds_idx_arr_dict_by_GeneSym(cds_info)

    GeneSym_list = []
    for ori_pos in pos_arr:
        pos = ori_pos + ADJ_REF_IDX

        for cds_arr in cds_info:
            GeneSym = cds_arr[0]
            # if GeneSym in GeneSym_list:
            #     continue
            start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)
            for i in range(len(start_idx_arr)):
                if start_idx_arr[i] <= pos <= end_idx_arr[i]:
                    GeneSym_list.append(GeneSym)
                    # break

    print(GeneSym_list)

    print(logic.make_complement_string('CCACACTGGGGCTCCCACTACTGCGAGGAGTGACCCACGAAGGCCACAGAGATGGCCGGGGCTTCGGTGAAGGTGGCGGTGCGGGTCCGCCCCTTCAATTCCCGGGAAATGAGCCGTGACTC')[::-1])
    print('CCACACTGGGGCTCCCACTACTGCGAGGAGTGACCCACGAAGGCCACAGAGATGGCCGGGGCTTCGGTGAAGGTGGCGGTGCGGGTCCGCCCCTTCAATTCCCGGGAAATGAGCCGTGACTC'[::-1])



if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    multi_processing_ClinVar_by_all_cds()
    # test()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
