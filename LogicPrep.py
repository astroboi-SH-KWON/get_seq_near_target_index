import re

class LogicPreps:

    def __init__(self):
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"

    def sort_list_by_ele(self, data_list, ele_idx, up_down_flag=True):
        result_list = []
        for tmp_arr in sorted(data_list, key=lambda tmp_arr: tmp_arr[ele_idx], reverse=up_down_flag):
            result_list.append(tmp_arr)
        return result_list

    """
    cds_file = [['', '', 'num_of_chromosome', '-', '14361', '29370', '29370', '29370', 'num_of_exon', 'start_idx_list', 'end_idx_list']
                , ['WASH7P', 'NR_024540', 'chr1', '-', '14361', '29370', '29370', '29370', '11', '14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320,', '14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370,']
                ]
    """
    def get_whole_idx_num_by_chr(self, cds_file):
        result_dict = {}
        for cds_info_arr in cds_file:
            chr_key = cds_info_arr[2]
            strand = cds_info_arr[3]
            start_dix_arr = re.findall('\d+', cds_info_arr[-2])
            end_dix_arr = re.findall('\d+', cds_info_arr[-1])

            tmp_list = []
            for idx in range(len(start_dix_arr)):
                nxt_idx = int(start_dix_arr[idx])
                end_idx = int(end_dix_arr[idx])
                tmp_list.extend([tmp_idx for tmp_idx in range(nxt_idx, end_idx)])

            if chr_key in result_dict:
                result_dict[chr_key].append(tmp_list)
            else:
                result_dict.update({chr_key: [tmp_list]})

        return result_dict

    def complement_char(self, ch):
        return {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[ch]

    def make_complement_string(self, trgt_seq):
        comp_seq = ""
        for ch in trgt_seq:
            comp_seq += self.complement_char(ch)
        return comp_seq

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




