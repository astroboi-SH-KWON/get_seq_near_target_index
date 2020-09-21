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
MUT_INFO = "20200909_ClinVar_hg38.txt"

TRGT_IDX = 686651



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


def read_FASTA():
    for seq_record in SeqIO.parse(REF_DIR + "chr1.fa", "fasta"):
        anti_seq = seq_record.seq.complement()
        print(seq_record.seq[TRGT_IDX - 1:TRGT_IDX - 1 + len("ATTTG")], " -1")
        print(anti_seq[TRGT_IDX - 1:TRGT_IDX - 1 + len("ATTTG")], " -1 complement")
        print(seq_record.seq[TRGT_IDX:TRGT_IDX + len("ATTTG")])
        print(anti_seq[TRGT_IDX:TRGT_IDX + len("ATTTG")])

def test():
    tmp_str = "PERM1	NM_001291366	chr1	-	975198	982117	976171	981029	4	975198,976498,978880,982064,	976269,976624,981047,982117,\n"
    tmp_str += "PERM1	NM_001291367	chr1	-	975198	982117	976171	982117	5	975198,976498,978880,981136,982064,	976269,976624,980657,981173,982117,\n"
    tmp_str += "SAMD11	NM_152486	chr1	+	925730	944574	925941	944153	14	925730,925921,930154,931038,935771,939039,939274,941143,942135,942409,942558,943252,943697,943907,	925800,926013,930336,931089,935896,939129,939460,941306,942251,942488,943058,943377,943808,944574,\n"
    logic_prep = LogicPrep.LogicPreps()

    cds_list = []
    for tmp_arr in tmp_str.split('\n')[:-1]:
        cds_list.append(tmp_arr.split('\t'))

    for cds_arr in cds_list:
        print(cds_arr)
        start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)
        print(start_idx_arr, " : start_idx_arr")
        print(end_idx_arr, " : end_idx_arr")
        print([end_idx_arr[i] - start_idx_arr[i] for i in range(len(start_idx_arr))], " : end_idx_arr[i] - start_idx_arr[i]")

def test1():
    logic = Logic.Logics()
    logic.exist_another_orf_end_codon_in_cds_seq("ABCD", False)







if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # main()
    # read_FASTA()
    test1()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
