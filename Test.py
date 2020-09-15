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

TRGT_IDX = 26696421



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


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    # main()
    read_FASTA()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
