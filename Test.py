import time
import os
from Bio import SeqIO
from Bio import motifs
import gzip
from Bio import SwissProt
from Bio import ExPASy
import glob

import Util
import Logic
import LogicPrep
############### start to set env ################
WORK_DIR = os.getcwd() + "/"
# REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/dat/gzip/"
# REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/dat/"
REF_DIR = "D:/000_WORK/000_reference_path/human/hg38/Splited/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
CDS_INFO = "hg38_refFlat_full.txt"
# MUT_INFO = "20200909_ClinVar_hg38.txt"
MUT_INFO = "200907_Dominant filter.txt"

TRGT_IDX = 686651

WIN_SIZE = [60, 60]
INIT_BY_PAM = [
    ['SaCas9', 'NNGRRT', 43, 3, WIN_SIZE]
    , ['SaCas9_KKH', 'NNNRRT', 43, 3, WIN_SIZE]
    , ['SaCas9_NNG', 'NNG', 43, 3, WIN_SIZE]
    , ['St1Cas9', 'NNRGAA', 41, 3, WIN_SIZE]
    , ['Nm1Cas9', 'NNNNGATT', 45, 3, WIN_SIZE]
    , ['Nm2Cas9', 'NNNNCC', 44, 3, WIN_SIZE]
    , ['CjCas9', 'NNNNRYAC', 44, 3, WIN_SIZE]
]

############### end setting env #################



def main():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()

    # CHROM	POS	ID	REF	ALT	mut_length	CLNVC	CLNSIG
    # [['1', '930188', '846933', 'G', 'A', '1', 'substitution', 'Uncertain_significance'],...]
    mut_list = util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t")
    # seq_idx - 0
    cds_file = util.read_csv_ignore_N_line(WORK_DIR + "input/" + CDS_INFO, "\t", 0)
    for idx in range(len(cds_file[:4])):
        print(cds_file[idx])
    cds_dict = logic_prep.get_whole_idx_num_by_chr(cds_file[:4])

    print(cds_dict)


def test_SY():
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics()
    util = Util.Utils()

    mut_list = util.read_csv_ignore_N_line(WORK_DIR + "input/" + MUT_INFO, "\t")[20000:]

    win_arr = WIN_SIZE
    for mut_arr in mut_list:
        tmp_arr = []
        tmp_arr.extend(mut_arr)
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

        for init in INIT_BY_PAM:
            logic.get_matched_pam_clvg_p_seq_dict(p_seq, pos, win_arr, ref_p_seq, init[1], init[2], init[3])

def test():
    tmp_str = "OR4F5	NM_001005484	chr1	+	69090	70008	69090	70008	1	69090,	70008,\n"
    tmp_str += "PERM1	NM_001291366	chr1	-	975198	982117	976171	981029	4	975198,976498,978880,982064,	976269,976624,981047,982117,\n"
    tmp_str += "PERM1	NM_001291367	chr1	-	975198	982117	976171	982117	5	975198,976498,978880,981136,982064,	976269,976624,980657,981173,982117,\n"
    tmp_str += "SAMD11	NM_152486	chr1	+	925730	944574	925941	944153	14	925730,925921,930154,931038,935771,939039,939274,941143,942135,942409,942558,943252,943697,943907,	925800,926013,930336,931089,935896,939129,939460,941306,942251,942488,943058,943377,943808,944574,\n"
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    cds_list = []
    for tmp_arr in tmp_str.split('\n')[:-1]:
        cds_list.append(tmp_arr.split('\t'))

    for cds_arr in cds_list:
        print(cds_arr)
        start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)
        idx_list = logic_prep.get_idx_num_frm_strt_to_end_list(start_idx_arr, end_idx_arr)
        p_seq, m_seq = util.read_file_by_biopython(REF_DIR + cds_arr[2] + ".fa", "fasta")
        print(p_seq[int(cds_arr[6]):int(cds_arr[7])])
        # p_cds_seq = logic_prep.get_seq_by_idx_arr(p_seq, start_idx_arr, end_idx_arr)
        # tmp_idx = ""
        # tmp_str = ""
        # for idx in idx_list:
        #     tmp_idx += "0"
        #     tmp_str += p_seq[idx]
        # print(len(tmp_idx), "len(tmp_idx)")
        # print(len(tmp_str), "len(tmp_str)")
        # print(tmp_str, "tmp_str")
        # print(p_cds_seq, "p_cds_seq")
        # if tmp_str == p_cds_seq:
        #     print("True")
        # else:
        #     print("False")

def test2():
    for i in range(5 - 3):
        print(i)




if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # main()
    test2()
    # test_SY()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
