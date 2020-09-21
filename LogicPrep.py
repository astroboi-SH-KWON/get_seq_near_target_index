import re

class LogicPreps:

    def __init__(self):
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"
        self.strt_NM_ = "NM_"

    def sort_list_by_ele(self, data_list, ele_idx, up_down_flag=True):
        result_list = []
        for tmp_arr in sorted(data_list, key=lambda tmp_arr: tmp_arr[ele_idx], reverse=up_down_flag):
            result_list.append(tmp_arr)
        return result_list

    """
    cds_file = [['GeneSym', 'NMID', 'Chrom', 'Strand', 'Transcript_Start', 'End', 'ORFStart', 'End', '#exon', 'start_idx_arr', 'end_idx_arr']
                , ['WASH7P', 'NR_024540', 'chr1', '-', '14361', '29370', '29370', '29370', '11', '14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320,', '14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370,']
                ]
    """
    def get_whole_idx_num_by_chr(self, cds_file):
        result_dict = {}
        for cds_info_arr in cds_file:
            chr_key = cds_info_arr[2]
            strand = cds_info_arr[3]
            print(strand)

            start_idx_arr, end_idx_arr = self.get_orf_strt_end_idx_arr(cds_info_arr)

            tmp_list = self.get_idx_num_frm_strt_to_end_list(start_idx_arr, end_idx_arr)

            if chr_key in result_dict:
                result_dict[chr_key].append(tmp_list)
            else:
                result_dict.update({chr_key: [tmp_list]})

        return result_dict

    def get_orf_strt_end_idx_arr(self, cds_arr):
        orf_strt_idx = int(cds_arr[6])
        orf_end_idx = int(cds_arr[7])

        start_idx_arr = [orf_strt_idx]
        start_idx_arr.extend([int(cds_idx) for cds_idx in re.findall('\d+', cds_arr[9])[1:]])
        end_idx_arr = [int(cds_idx) for cds_idx in re.findall('\d+', cds_arr[10])[:-1]]
        end_idx_arr.append(orf_end_idx)

        # check UTR in first and last codon
        sorted_start_idx_arr = sorted(start_idx_arr)
        sorted_end_idx_arr = sorted(end_idx_arr)
        strt_idx = sorted_start_idx_arr.index(orf_strt_idx)
        end_idx = sorted_end_idx_arr.index(orf_end_idx)

        return sorted_start_idx_arr[strt_idx: end_idx + 1], sorted_end_idx_arr[strt_idx: end_idx + 1]

    def get_idx_num_frm_strt_to_end_list(self, start_idx_arr, end_idx_arr):
        result_list = []
        for idx in range(len(start_idx_arr)):
            nxt_idx = int(start_idx_arr[idx])
            end_idx = int(end_idx_arr[idx])
            result_list.extend([tmp_idx for tmp_idx in range(nxt_idx, end_idx)])
        return result_list

    def get_seq_by_idx_arr(self, full_seq, start_idx_arr, end_idx_arr):
        result_seq = ""
        for idx in range(len(start_idx_arr)):
            nxt_idx = int(start_idx_arr[idx])
            end_idx = int(end_idx_arr[idx])
            result_seq = full_seq[nxt_idx: end_idx]

        return result_seq


    def add_result_seq_to_arr(self, mut_arr, data_dict):
        tmp_str = ""
        for tmp_arr in data_dict.values():
            for tmp_seq in tmp_arr:
                tmp_str += tmp_seq + " , "
        mut_arr.append(tmp_str)

    def merge_multi_list(self, pool_list):
        result_list = []
        for split_list in pool_list:
            result_list.extend(split_list)

        return result_list

    def filter_out_NON_NM_id_in_cds_list(self, cds_list):
        result_list = []
        gene_sym_dict = {}

        file_nm_arr = ['chrX', 'chrY']
        for f_num in range(1, 23):
            file_nm_arr.append("chr" + str(f_num))

        for cds_arr in cds_list:
            gene_sym = cds_arr[0]
            trgt_nm_id = cds_arr[1]

            file_nm = cds_arr[2]

            # filter out if 'NM_' in trgt_nm_id
            if not trgt_nm_id.startswith(self.strt_NM_):
                continue
            # filter out if chromosome file isn't verified like 'chr19_KI270917v1_alt'...
            if file_nm not in file_nm_arr:
                continue

            # make nm_id numerical for sorting
            nm_num = int(cds_arr[1].lstrip(self.strt_NM_).lstrip("0"))
            cds_arr.append(nm_num)
            if gene_sym in gene_sym_dict:
                gene_sym_dict[gene_sym].append(cds_arr)
            else:
                gene_sym_dict.update({gene_sym: [cds_arr]})

        # get the lowest num of nm_id in each gene_sym
        for cds_list in gene_sym_dict.values():
            sorted_cds_list = self.sort_list_by_ele(cds_list, -1, False)
            result_list.append(sorted_cds_list[0])

        return result_list




