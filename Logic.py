from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62

import LogicPrep
class Logics:
    def __init__(self):
        self.strt_cd_arr = ['ATG']
        self.end_cd_arr = ['TGA', 'TAG', 'TAA']

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

    def get_matched_pam_p_seq_dict(self, p_seq, pos, win_arr, ref_seq, pam, len_f_pam, len_b_pam):
        result_dict = {}
        len_pam = len(pam)
        len_ref_seq = len(ref_seq)

        p_ori_win_f_seq = p_seq[pos - win_arr[0]: pos]
        p_ori_win_b_seq = p_seq[pos + len(ref_seq): pos + len(ref_seq) + win_arr[1]]
        p_ori_win_seq = p_ori_win_f_seq + ref_seq + p_ori_win_b_seq

        p_first_pam_pos = len(p_ori_win_f_seq) - len_pam + 1

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

    def get_matched_pam_m_seq_dict(self, m_seq, pos, win_arr, ref_seq, pam, len_f_pam, len_b_pam):
        result_dict = {}
        len_pam = len(pam)
        len_ref_seq = len(ref_seq)

        m_ori_win_f_seq = m_seq[pos - win_arr[0]: pos]
        m_ori_win_b_seq = m_seq[pos + len(ref_seq): pos + len(ref_seq) + win_arr[1]]
        m_ori_win_seq = m_ori_win_f_seq + ref_seq + m_ori_win_b_seq

        m_first_pam_pos = len(m_ori_win_f_seq) - len_pam + 1

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

        for mut_arr in mut_list:
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
