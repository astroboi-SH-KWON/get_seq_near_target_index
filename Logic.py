from Bio import SeqIO

import LogicPrep
import Util
class Logics:
    def __init__(self):
        self.strt_cd_arr = ['ATG']
        self.end_cd_arr = ['TGA', 'TAG', 'TAA']
        self.len_clvg = 3
        self.len_spcr = 21
        self.len_contxt = 22

    def complement_char(self, ch):
        complement_char_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        try:
            return complement_char_dict[ch]
        except:
            print("complement_char : [" + ch + "]")
            raise Exception

    def make_complement_string(self, trgt_seq):
        comp_seq = ""
        for ch in trgt_seq:
            try:
                comp_seq += self.complement_char(ch)
            except:
                print(trgt_seq, 'make_complement_string')
                raise Exception
        return comp_seq

    """
    checkSeqByChar : match sequences by char
    :param
        seq_char :
        target_char : 
    :return
        boolean
    """
    def checkSeqByChar(self, seq_char, target_char):
        flag = False
        if target_char == 'N':
            return True
        elif target_char in 'ACGTU':
            if seq_char == target_char:
                return True
        elif target_char == 'R':
            if seq_char in 'AG':
                return True
        elif target_char == 'Y':
            if seq_char in 'CT':
                return True
        """
        add more rules of "ACGTU"
        """

        return flag

    def checkSeqByChar_SY(self, seq_char, target_char):
        flag = False
        if target_char == 'N' or seq_char == 'N':
            return True
        elif target_char in 'ACGTU':
            if seq_char == target_char:
                return True

        elif target_char == 'M':
            if seq_char in 'AC':
                return True
        elif target_char == 'R':
            if seq_char in 'AG':
                return True
        elif target_char == 'W':
            if seq_char in 'AT':
                return True
        elif target_char == 'S':
            if seq_char in 'CG':
                return True
        elif target_char == 'Y':
            if seq_char in 'CT':
                return True
        elif target_char == 'K':
            if seq_char in 'GT':
                return True
        elif target_char == 'V':
            if seq_char != 'T':
                return True
        elif target_char == 'H':
            if seq_char != 'G':
                return True
        elif target_char == 'D':
            if seq_char != 'C':
                return True
        elif target_char == 'B':
            if seq_char != 'A':
                return True
        """
        add more rules of "ACGTU"
        """

        return flag

    """
    match : match sequence with "same length" strings
    :param
        i : index of seq
        seq_str : targeted DNA/RNA sequence 
        rule_str : rules with "ACGTU", "N", "R",...
    :return
        boolean
    """
    def match(self, i, seq_str, rule_str):
        if len(seq_str) == i:
            return True
        if self.checkSeqByChar(seq_str[i], rule_str[i]):
            return self.match(i + 1, seq_str, rule_str)
        else:
            return False

    def match_SY(self, i, seq_str, rule_str):
        if len(seq_str) == i:
            return True
        if self.checkSeqByChar_SY(seq_str[i], rule_str[i]):
            return self.match_SY(i + 1, seq_str, rule_str)
        else:
            return False

    def get_matched_pam_p_seq_dict(self, p_seq, pos, win_arr, pos_seq, pam, len_f_pam, len_b_pam):
        result_dict = {}
        len_pam = len(pam)
        len_ref_seq = len(pos_seq)

        p_ori_win_f_seq = p_seq[pos - win_arr[0]: pos]
        p_ori_win_b_seq = p_seq[pos + len(pos_seq): pos + len(pos_seq) + win_arr[1]]
        p_ori_win_seq = p_ori_win_f_seq + pos_seq + p_ori_win_b_seq

        p_first_pam_pos = len(p_ori_win_f_seq) - len_pam + 1

        # check PAM in only pos_seq part
        for i in range(len_pam + len_ref_seq - 1):
            p_ref_pam_seq = p_ori_win_seq[i + p_first_pam_pos: i + p_first_pam_pos + len_pam]
            if self.match_SY(0, p_ref_pam_seq, pam):
                p_seq_f_pam = p_ori_win_seq[i + p_first_pam_pos - len_f_pam: i + p_first_pam_pos]
                p_seq_b_pam = p_ori_win_seq[i + p_first_pam_pos + len_pam: i + p_first_pam_pos + len_pam + len_b_pam]
                if p_seq_f_pam in result_dict:
                    result_dict[p_seq_f_pam].append(p_seq_f_pam + p_ref_pam_seq + p_seq_b_pam)
                else:
                    result_dict.update({p_seq_f_pam: [p_seq_f_pam + p_ref_pam_seq + p_seq_b_pam]})

        return result_dict, p_ori_win_seq

    def get_matched_pam_m_seq_dict(self, m_seq, pos, win_arr, pos_seq, pam, len_f_pam, len_b_pam):
        result_dict = {}
        len_pam = len(pam)
        len_ref_seq = len(pos_seq)

        m_ori_win_f_seq = m_seq[pos - win_arr[0]: pos]
        m_ori_win_b_seq = m_seq[pos + len(pos_seq): pos + len(pos_seq) + win_arr[1]]
        m_ori_win_seq = m_ori_win_f_seq + pos_seq + m_ori_win_b_seq

        m_first_pam_pos = len(m_ori_win_f_seq) - len_pam + 1

        # check PAM in only pos_seq part
        for i in range(len_pam + len_ref_seq - 1):
            m_ref_pam_seq = m_ori_win_seq[i + m_first_pam_pos: i + m_first_pam_pos + len_pam]
            if self.match_SY(0, m_ref_pam_seq, pam[::-1]):
                m_seq_f_pam = m_ori_win_seq[i + m_first_pam_pos + len_pam: i + m_first_pam_pos + len_pam + len_f_pam]
                m_seq_b_pam = m_ori_win_seq[i + m_first_pam_pos - len_b_pam: i + m_first_pam_pos]
                if m_seq_f_pam[::-1] in result_dict:
                    result_dict[m_seq_f_pam[::-1]].append((m_seq_b_pam + m_ref_pam_seq + m_seq_f_pam)[::-1])
                else:
                    result_dict.update({m_seq_f_pam[::-1]: [(m_seq_b_pam + m_ref_pam_seq + m_seq_f_pam)[::-1]]})

        return result_dict, m_ori_win_seq[::-1]

    def get_matched_pam_clvg_p_seq_dict(self, p_seq, pos, win_arr, pos_seq, pam, len_f_pam, len_b_pam):
        result_dict = {}
        len_pam = len(pam)
        len_ref_seq = len(pos_seq)

        p_ori_win_f_seq = p_seq[pos - win_arr[0]: pos]
        p_ori_win_b_seq = p_seq[pos + len(pos_seq): pos + len(pos_seq) + win_arr[1]]
        p_ori_win_seq = p_ori_win_f_seq + pos_seq + p_ori_win_b_seq

        p_first_pam_pos = len(p_ori_win_f_seq) - len_pam + 1

        # check PAM in only pos_seq part
        for i in range(len_pam + len_ref_seq - 1):
            p_ref_pam_seq = p_ori_win_seq[i + p_first_pam_pos: i + p_first_pam_pos + len_pam]
            if self.match_SY(0, p_ref_pam_seq, pam):
                p_seq_f_pam = p_ori_win_seq[i + p_first_pam_pos - len_f_pam: i + p_first_pam_pos]
                p_seq_b_pam = p_ori_win_seq[i + p_first_pam_pos + len_pam: i + p_first_pam_pos + len_pam + len_b_pam]

                pos_clvg = i + pos - len_pam + 1 - self.len_clvg
                if p_seq_f_pam in result_dict:
                    result_dict[p_seq_f_pam].append([p_seq_f_pam + p_ref_pam_seq + p_seq_b_pam, pos_clvg])
                else:
                    result_dict.update({p_seq_f_pam: [[p_seq_f_pam + p_ref_pam_seq + p_seq_b_pam, pos_clvg]]})

        return result_dict, p_ori_win_seq

    def get_matched_pam_clvg_m_seq_dict(self, m_seq, pos, win_arr, pos_seq, pam, len_f_pam, len_b_pam):
        result_dict = {}
        len_pam = len(pam)
        len_ref_seq = len(pos_seq)

        m_ori_win_f_seq = m_seq[pos - win_arr[0]: pos]
        m_ori_win_b_seq = m_seq[pos + len(pos_seq): pos + len(pos_seq) + win_arr[1]]
        m_ori_win_seq = m_ori_win_f_seq + pos_seq + m_ori_win_b_seq

        m_first_pam_pos = len(m_ori_win_f_seq) - len_pam + 1

        # check PAM in only pos_seq part
        for i in range(len_pam + len_ref_seq - 1):
            m_ref_pam_seq = m_ori_win_seq[i + m_first_pam_pos: i + m_first_pam_pos + len_pam]
            if self.match_SY(0, m_ref_pam_seq, pam[::-1]):
                m_seq_f_pam = m_ori_win_seq[i + m_first_pam_pos + len_pam: i + m_first_pam_pos + len_pam + len_f_pam]
                m_seq_b_pam = m_ori_win_seq[i + m_first_pam_pos - len_b_pam: i + m_first_pam_pos]

                pos_clvg = i + pos + 1 + self.len_clvg
                if m_seq_f_pam[::-1] in result_dict:
                    result_dict[m_seq_f_pam[::-1]].append([(m_seq_b_pam + m_ref_pam_seq + m_seq_f_pam)[::-1], pos_clvg])
                else:
                    result_dict.update({m_seq_f_pam[::-1]: [[(m_seq_b_pam + m_ref_pam_seq + m_seq_f_pam)[::-1], pos_clvg]]})

        return result_dict, m_ori_win_seq[::-1]

    def get_matched_pam_p_seq_list(self, result_list, p_seq, pos, pam, win_pam, seq_win_size, chr):
        win_pam_f_pos = win_pam[0]
        win_pam_b_pos = win_pam[1]

        len_f_pam = seq_win_size[0]
        len_b_pam = seq_win_size[1]

        p_seq_pam_win = p_seq[pos - win_pam_f_pos: pos + win_pam_b_pos]
        p_pam_pos = pos - win_pam_f_pos

        for i in range(len(p_seq_pam_win) - len(pam) + 1):
            candi_pam = p_seq_pam_win[i: i + len(pam)]
            if self.match(0, candi_pam, pam):
                result_list.append(
                    [chr, candi_pam, p_seq[p_pam_pos - len_f_pam: p_pam_pos + len(pam) + len_b_pam], p_pam_pos, '+'])
            p_pam_pos += 1

    def get_matched_pam_m_seq_list(self, result_list, m_seq, pos, pam, win_pam, seq_win_size, chr):
        win_pam_f_pos = win_pam[0]
        win_pam_b_pos = win_pam[1]

        len_f_pam = seq_win_size[0]
        len_b_pam = seq_win_size[1]

        m_seq_pam_win = m_seq[pos - win_pam_b_pos: pos + win_pam_f_pos]
        m_pam_pos = pos - win_pam_b_pos

        for i in range(len(m_seq_pam_win) - len(pam) + 1):
            candi_pam = m_seq_pam_win[i: i + len(pam)]
            if self.match(0, candi_pam, pam[::-1]):
                result_list.append(
                    [chr, candi_pam[::-1], m_seq[m_pam_pos - len_b_pam: m_pam_pos + len(pam) + len_f_pam][::-1],
                     m_pam_pos + len(pam), '-'])
            m_pam_pos += 1

    def remove_dict_val_by_key(self, trgt_dict, keys):
        for tmp_key in keys:
            if tmp_key in trgt_dict:
                trgt_dict.pop(tmp_key)

    def get_seq_by_pam_after_mut(self, path, mut_list, win_arr, init):
        logic_prep = LogicPrep.LogicPreps()

        pam_arr = init[1]
        len_f_pam_arr = init[2]
        len_b_pam_arr = init[3]
        adj_ref_idx = init[4]

        for mut_arr in mut_list:
            chr_num = mut_arr[0]
            pos = int(mut_arr[1]) + adj_ref_idx
            ref_p_seq = mut_arr[3]
            alt_p_seq = mut_arr[4]

            ref_m_seq = ""
            alt_m_seq = ""
            try:
                ref_m_seq += self.make_complement_string(ref_p_seq)
                if alt_p_seq == '.':
                    alt_p_seq = ""
                else:
                    alt_m_seq += self.make_complement_string(alt_p_seq)
            except Exception as err:
                print("make_complement_string ::: ", err)
                print(ref_p_seq, " : ref_p_seq")
                print(alt_p_seq, " : alt_p_seq")
                print(str(mut_arr))

            seq_record = SeqIO.read(path + "chr" + chr_num + ".fa", "fasta")
            p_seq = str(seq_record.seq).upper()
            m_seq = str(seq_record.seq.complement()).upper()

            ori_win_flag = True

            for idx in range(len(pam_arr)):
                pam = pam_arr[idx]
                len_f_pam = len_f_pam_arr[idx]
                len_b_pam = len_b_pam_arr[idx]

                ref_p_dict, p_ori_win_seq = self.get_matched_pam_p_seq_dict(p_seq, pos, win_arr, ref_p_seq, pam,
                                                                             len_f_pam, len_b_pam)
                ref_m_dict, m_ori_win_seq = self.get_matched_pam_m_seq_dict(m_seq, pos, win_arr, ref_m_seq, pam,
                                                                             len_f_pam, len_b_pam)

                mut_p_dict, _ = self.get_matched_pam_p_seq_dict(p_seq, pos, win_arr, alt_p_seq, pam, len_f_pam,
                                                                 len_b_pam)
                mut_m_dict, _ = self.get_matched_pam_m_seq_dict(m_seq, pos, win_arr, alt_m_seq, pam, len_f_pam,
                                                                 len_b_pam)

                self.remove_dict_val_by_key(mut_p_dict, ref_p_dict.keys())
                self.remove_dict_val_by_key(mut_m_dict, ref_m_dict.keys())

                if ori_win_flag:
                    mut_arr.append(p_ori_win_seq + " , " + m_ori_win_seq)
                    ori_win_flag = False

                logic_prep.add_result_seq_to_arr(mut_arr, mut_p_dict)
                logic_prep.add_result_seq_to_arr(mut_arr, mut_m_dict)

    def is_all_bigger_than_cut_off(self, int_list, cut_off=0):
        return all(i > cut_off for i in int_list)

    def exist_another_orf_end_codon_in_cds_seq(self, cds_seq, flag=True):
        if not flag:
            cds_seq = cds_seq[::-1]

        codons_wo_last_cd_list = [cds_seq[i * 3:(i + 1) * 3] for i in range(int(len(cds_seq) / 3 - 1))]
        for end_cd in self.end_cd_arr:
            if end_cd in codons_wo_last_cd_list:
                return True
        return False

    def get_len_each_codon_list(self, start_idx_arr, end_idx_arr):
        return [end_idx_arr[i] - start_idx_arr[i] for i in range(len(start_idx_arr))]

    def get_seq_pam_in_orf(self, key, cds_dict, seq_idx_dict, result_list):
        logic_prep = LogicPrep.LogicPreps()

        print("key : ", key, " , start >>> get_seq_pam_in_orf")
        cds_list = cds_dict[key]
        seq_idx_list = seq_idx_dict[key]

        for cds_arr in cds_list:
            gene_sym = cds_arr[0]
            nm_id = cds_arr[1]
            chr_nm = cds_arr[2]
            strand = cds_arr[3]
            orf_strt_pos = int(cds_arr[6])
            orf_end_pos = int(cds_arr[7])

            start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)
            idx_list = logic_prep.get_idx_num_frm_strt_to_end_list(start_idx_arr, end_idx_arr)

            for seq_idx_arr in seq_idx_list:
                pam_full_idx = int(seq_idx_arr[3])
                if pam_full_idx in idx_list:
                    pam_idx = idx_list.index(pam_full_idx)

                    if strand != seq_idx_arr[4]:
                        continue
                    tmp_arr = []
                    tmp_arr.extend(seq_idx_arr)
                    if strand == '+':
                        orf_missing = pam_idx % 3
                        tmp_arr.extend([orf_missing, gene_sym, nm_id])
                    else:
                        orf_missing = (len(idx_list) - pam_idx - 1) % 3
                        tmp_arr.extend([orf_missing, gene_sym, nm_id])
                    result_list.append(tmp_arr)

        print("DONE : ", key, " >>> get_seq_pam_in_orf")

    def get_pos_ratio_in_cds(self, strand, idx, idx_list, pam_seq):
        if strand == '+':
            return (idx + 1 - self.len_clvg) / len(idx_list)
        else:
            return (len(idx_list) - idx - len(pam_seq) - self.len_clvg) / len(idx_list)

    def get_pos_ratio_in_cds_by_full_idx(self, strand, cdn_idx, idx_list, pam_seq):
        idx = idx_list.index(cdn_idx)
        if strand == '+':
            return (idx + 1 - self.len_clvg) / len(idx_list)
        else:
            return (len(idx_list) - idx - len(pam_seq) - self.len_clvg) / len(idx_list)

    def get_seq_clvg_pos_in_cds_w_whole_gene(self, ref_dir, key, init, result_list):
        def_nm = "get_seq_clvg_pos_in_cds_w_whole_gene"
        print("key : ", key, " , start >>> ", def_nm)
        pam_seq = init[1]
        len_f_pam = init[2]
        len_b_pam = init[3]
        win_arr = init[4]
        ratio_f = init[5][0]
        ratio_b = init[5][1]
        cds_list = init[6][key]

        logic_prep = LogicPrep.LogicPreps()
        util = Util.Utils()
        p_seq, m_seq = util.read_file_by_biopython(ref_dir + key + ".fa", "fasta")

        for cds_arr in cds_list:
            gene_sym = cds_arr[0]
            nm_id = cds_arr[1]
            chr_nm = cds_arr[2]
            strand = cds_arr[3]

            start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)
            idx_list = logic_prep.get_idx_num_frm_strt_to_end_list(start_idx_arr, end_idx_arr)

            for clvg_site in idx_list:
                # filter out the edge index of codons
                if self.is_in_list(clvg_site, start_idx_arr) and self.is_in_list(clvg_site, end_idx_arr):
                    continue

                pos_ratio_cds = self.get_pos_ratio_in_cds_by_full_idx(strand, clvg_site, idx_list, pam_seq)
                if ratio_f < pos_ratio_cds < ratio_b:
                    trgt_p_pam = p_seq[clvg_site + self.len_clvg: clvg_site + self.len_clvg + len(pam_seq)]
                    trgt_m_pam = m_seq[clvg_site - self.len_clvg - len(pam_seq): clvg_site - self.len_clvg][::-1]

                    # check + strand
                    if self.match_SY(0, trgt_p_pam, pam_seq):
                        trgt_full_p_seq = p_seq[clvg_site + self.len_clvg - len_f_pam: clvg_site + self.len_clvg + len(
                            pam_seq) + len_b_pam]
                        trgt_p_seq = p_seq[
                                     clvg_site + self.len_clvg - self.len_spcr - self.len_contxt: clvg_site + self.len_clvg + len(
                                         pam_seq) + len_b_pam]
                        result_list.append(
                            [chr_nm, gene_sym, nm_id, strand, clvg_site + 1, '+', trgt_p_seq[:self.len_contxt],
                             trgt_p_seq[self.len_contxt: self.len_contxt + self.len_spcr], trgt_p_pam,
                             trgt_p_seq[- len_b_pam:], trgt_full_p_seq, pos_ratio_cds])

                    # check - strand
                    if self.match_SY(0, trgt_m_pam, pam_seq):
                        trgt_full_m_seq = m_seq[clvg_site - self.len_clvg - len(
                            pam_seq) - len_b_pam: clvg_site - self.len_clvg + len_f_pam][::-1]
                        trgt_m_seq = m_seq[clvg_site - self.len_clvg - len(
                            pam_seq) - len_b_pam: clvg_site - self.len_clvg + self.len_spcr + self.len_contxt][::-1]
                        result_list.append(
                            [chr_nm, gene_sym, nm_id, strand, clvg_site + 1, '-', trgt_m_seq[:self.len_contxt],
                             trgt_m_seq[self.len_contxt: self.len_contxt + self.len_spcr], trgt_p_pam,
                             trgt_m_seq[- len_b_pam:], trgt_full_m_seq, pos_ratio_cds])

        print("DONE : ", key, " >>> ", def_nm)

    def get_seq_clvg_pos_in_cds(self, ref_dir, key, init, result_list):
        def_nm = "get_seq_clvg_pos_in_cds"
        print("key : ", key, " , start >>> ", def_nm)
        otholog = init[0]
        pam_seq = init[1]
        len_f_pam = init[2]
        len_b_pam = init[3]
        win_arr = init[4]
        ratio_f = init[5][0]
        ratio_b = init[5][1]
        cds_list = init[6][key]
        mut_list = init[7][key.replace("chr", "")]

        logic_prep = LogicPrep.LogicPreps()

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
                ref_m_seq += self.make_complement_string(ref_p_seq)
                if alt_p_seq == '.':
                    alt_p_seq = ""
                else:
                    alt_m_seq += self.make_complement_string(alt_p_seq)
            except Exception as err:
                print(def_nm + " ==> make_complement_string ::: ", err)
                print(ref_p_seq, " : ref_p_seq")
                print(alt_p_seq, " : alt_p_seq")
                print(str(mut_arr))

            seq_record = SeqIO.read(ref_dir + "chr" + chr_num + ".fa", "fasta")
            p_seq = str(seq_record.seq).upper()
            m_seq = str(seq_record.seq.complement()).upper()

            ref_p_dict, p_ori_win_seq = self.get_matched_pam_clvg_p_seq_dict(p_seq, pos, win_arr, ref_p_seq, pam_seq,
                                                                         len_f_pam, len_b_pam)
            ref_m_dict, m_ori_win_seq = self.get_matched_pam_clvg_m_seq_dict(m_seq, pos, win_arr, ref_m_seq, pam_seq,
                                                                         len_f_pam, len_b_pam)

            mut_p_dict, _ = self.get_matched_pam_clvg_p_seq_dict(p_seq, pos, win_arr, alt_p_seq, pam_seq, len_f_pam,
                                                             len_b_pam)
            mut_m_dict, _ = self.get_matched_pam_clvg_m_seq_dict(m_seq, pos, win_arr, alt_m_seq, pam_seq, len_f_pam,
                                                             len_b_pam)

            self.remove_dict_val_by_key(mut_p_dict, ref_p_dict.keys())
            self.remove_dict_val_by_key(mut_m_dict, ref_m_dict.keys())

            for cds_arr in cds_list:
                gene_sym = cds_arr[0]
                nm_id = cds_arr[1]
                chr_nm = cds_arr[2]
                strand = cds_arr[3]

                start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)
                idx_list = logic_prep.get_idx_num_frm_strt_to_end_list(start_idx_arr, end_idx_arr)

                for val_list in mut_p_dict.values():
                    for val_arr in val_list:
                        res_seq = val_arr[0]
                        pos_clvg = val_arr[1]

                        if pos_clvg in idx_list:
                            if strand == '+':
                                idx_clvg_cds = idx_list.index(pos_clvg) + 1
                                pos_ratio_cds = idx_clvg_cds / len(idx_list)
                                if ratio_f < pos_ratio_cds < ratio_b:
                                    result_list.append(
                                        [chr_nm, gene_sym, nm_id, pos_clvg, strand, res_seq[:22], res_seq[22:len_f_pam],
                                         res_seq[len_f_pam: - len_b_pam], res_seq[- len_b_pam:], res_seq,
                                         pos_ratio_cds])

                            else:
                                idx_clvg_cds = idx_list.index(pos_clvg)
                                pos_ratio_cds = (len(idx_list) - idx_clvg_cds) / len(idx_list)
                                if ratio_f < pos_ratio_cds < ratio_b:
                                    result_list.append(
                                        [chr_nm, gene_sym, nm_id, pos_clvg, strand, res_seq[:22], res_seq[22:len_f_pam],
                                         res_seq[len_f_pam: - len_b_pam], res_seq[- len_b_pam:], res_seq,
                                         pos_ratio_cds])

                for val_list in mut_m_dict.values():
                    for val_arr in val_list:
                        res_seq = val_arr[0]
                        pos_clvg = val_arr[1]

                        if pos_clvg in idx_list:
                            if strand == '+':
                                idx_clvg_cds = idx_list.index(pos_clvg) + 1
                                pos_ratio_cds = idx_clvg_cds / len(idx_list)
                                if ratio_f < pos_ratio_cds < ratio_b:
                                    result_list.append(
                                        [chr_nm, gene_sym, nm_id, pos_clvg, strand, res_seq[:22], res_seq[22:len_f_pam],
                                         res_seq[len_f_pam: - len_b_pam], res_seq[- len_b_pam:], res_seq,
                                         pos_ratio_cds])

                            else:
                                idx_clvg_cds = idx_list.index(pos_clvg)
                                pos_ratio_cds = (len(idx_list) - idx_clvg_cds) / len(idx_list)
                                if ratio_f < pos_ratio_cds < ratio_b:
                                    result_list.append(
                                        [chr_nm, gene_sym, nm_id, pos_clvg, strand, res_seq[:22], res_seq[22:len_f_pam],
                                         res_seq[len_f_pam: - len_b_pam], res_seq[- len_b_pam:], res_seq,
                                         pos_ratio_cds])

        print("DONE : ", key, " >>> ", def_nm)

    def get_num_exon(self, pos_clvg, start_idx_arr, end_idx_arr, flag=True):
        total_exon = len(start_idx_arr)
        for j in range(total_exon):
            strt_pos = start_idx_arr[j]
            end_pos = end_idx_arr[j]
            if strt_pos < pos_clvg < end_pos:
                if flag:
                    return j + 1
                else:
                    return total_exon - j

    def check_strand(self, p_seq, m_seq, trgt_seq):
        if p_seq == trgt_seq:
            return '+'
        elif m_seq == trgt_seq:
            return '-'
        else:
            return

    def is_in_list(self, idx, check_list):
        if idx in check_list:
            return True
        return False

    def check_seq_in_cds(self, cds_idx_list, trgt_pos):
        for cds_idx_arr in cds_idx_list:
            if trgt_pos in cds_idx_arr:
                return True
        return False

    # GeneSym	NMID	Chrom	Strand	Transcript_Start	End	ORFStart	End	#Exon	ExonS_list	ExonE_list
    # Xkr4	CCDS14803.1	chr1	-	3216021	3671347	3216021	3671347	3	3216021,3421701,3670551,	3216967,3421900,3671347,
    def get_cds_idx_arr_dict_by_GeneSym(self, cds_info):
        logic_prep = LogicPrep.LogicPreps()

        result_dict = {}
        for cds_arr in cds_info:
            gene_sym = cds_arr[0]
            start_idx_arr, end_idx_arr = logic_prep.get_orf_strt_end_idx_arr(cds_arr)
            idx_list = logic_prep.get_idx_num_frm_strt_to_end_list(start_idx_arr, end_idx_arr)
            if gene_sym in result_dict:
                result_dict[gene_sym].append(idx_list)
            else:
                result_dict.update({gene_sym: [idx_list]})
        return result_dict

    def is_on_cds(self, clvg_site, gene_list, cds_dict):
        for gen_nm in gene_list:
            cds_idx_list_arr = cds_dict[gen_nm]
            for cds_idx_list in cds_idx_list_arr:
                if clvg_site in cds_idx_list:
                    return True
        return False

    # STRAND, REF (2.3에서의 sequence), Guide context (22 nt + guide RNA + PAM + 3nt)
    def is_alt_guide_in_ref_seq(self, init_arr, val_arr):
        len_f_guide = init_arr[1]
        len_guide = init_arr[2]
        pam_rule_arr = init_arr[3]
        # len_b_pam = init_arr[4]

        ref_seq = val_arr[1]
        cntxt_seq = val_arr[2]

        for pam_seq in pam_rule_arr:
            len_pam = len(pam_seq)
            for i in range(len_guide, len(ref_seq) - len_pam):
                pam_fr_ref = ref_seq[i: i + len_pam]
                if self.match_SY(0, pam_fr_ref, pam_seq):
                    guide_fr_ref = ref_seq[i - len_guide: i]
                    guide_fr_cntxt = cntxt_seq[len_f_guide: len_f_guide + len_guide]
                    if guide_fr_cntxt == guide_fr_ref:
                        return True
        return False

