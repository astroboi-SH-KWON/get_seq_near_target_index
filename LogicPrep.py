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
            result_seq += full_seq[nxt_idx: end_idx]

        return result_seq


    def get_seq_arr_by_idx_arr(self, full_seq, start_idx_arr, end_idx_arr):
        result_list = []
        for idx in range(len(start_idx_arr)):
            nxt_idx = int(start_idx_arr[idx])
            end_idx = int(end_idx_arr[idx])
            result_list.append(full_seq[nxt_idx: end_idx])

        return result_list


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
            result_list.append(sorted_cds_list[0][:-1])

        return result_list

    def get_dict_from_list_by_ele_key(self, data_list, key_idx):
        result_dict = {}
        for val_arr in data_list:
            if val_arr[key_idx] in result_dict:
                result_dict[val_arr[key_idx]].append(val_arr)
            else:
                result_dict.update({val_arr[key_idx]: [val_arr]})
        return result_dict

    def get_data_with_trgt_strng(self, ccds_list, trgt_str, idx):
        return [cds_arr for cds_arr in ccds_list if cds_arr[idx] == trgt_str]

    def get_high_ccds_id_among_same_gen_id(self, ccds_list):
        gen_id_dict = {}
        result_list = []
        for ccds_arr in ccds_list:
            gene = ccds_arr[2]
            ccds_id = ''.join(x for x in ccds_arr[4] if x.isdigit())
            ccds_arr.append(ccds_id)
            if gene in gen_id_dict:
                gen_id_dict[gene].append(ccds_arr)
            else:
                gen_id_dict.update({gene: [ccds_arr]})

        for ccds_list in gen_id_dict.values():
            """
            sangyeon  오후 6:16 20201123
            그리고 CCDS id의 순선는 큰 숫자가 더 최근에 나온 data네요!
            """
            sorted_ccds_list = self.sort_list_by_ele(ccds_list, -1, True)
            result_list.append(sorted_ccds_list[0][:-1])
        return result_list


    """
    mouse_ccds form : 
        #chromosome	nc_accession	gene	gene_id	ccds_id	ccds_status	cds_strand	cds_from	cds_to	cds_locations	match_type
    hg38_refFlat form : 
        GeneSym   NMID    Chrom   Strand  Transcript_Start   End ORFStart    End #Exon   ExonS_list  ExonE_list
    
    GeneSyn = gene
    NMID = ccds_id (이 부분은 human을 filter하실 때, gene 하나당 NMID가 여러개였던걸 줄이신거면 동일한 것으로 생각됩니다)
    Chrom = #chromosome
    Strand = cds_strand
    Transcript_Start =
    End =
    ORF Start = cds_from
    End = cds_to
    #Exon, ExonS_list, ExonE_list = cds_locations
    """
    def transform_mouse_ccds_form_to_hg38_refFlat(self, ccds_list):
        result_list = []
        for ccds_arr in ccds_list:
            gene_sym = ccds_arr[2]
            nm_id = ccds_arr[4]
            chrom = 'chr' + str(ccds_arr[0])
            strand = ccds_arr[6]

            transcript_st = ccds_arr[7]
            transcript_en = ccds_arr[8]
            orf_st = ccds_arr[7]
            orf_en = ccds_arr[8]

            cds_locations_arr = re.findall('\d+', ccds_arr[9])
            num_exon = len(cds_locations_arr) // 2
            exon_s_list = ''
            for cds_loc in [cds_locations_arr[i] for i in range(len(cds_locations_arr)) if i % 2 == 0]:
                exon_s_list = exon_s_list + str(cds_loc) + ','
            exon_e_list = ''
            for cds_loc in [cds_locations_arr[i] for i in range(len(cds_locations_arr)) if i % 2 != 0]:
                exon_e_list = exon_e_list + str(cds_loc) + ','

            result_list.append(
                [gene_sym, nm_id, chrom, strand, transcript_st, transcript_en, orf_st, orf_en, num_exon, exon_s_list,
                 exon_e_list])

        return result_list
