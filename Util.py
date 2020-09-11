import glob
from Bio import SeqIO
import openpyxl
import os

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    def get_seq_record_from_genbank(self, path):
        return SeqIO.read(path, "genbank")

    def make_row(self, sheet, row, data_arr, col=1):
        for idx in range(len(data_arr)):
            sheet.cell(row=row, column=(col + idx), value=data_arr[idx])

    def make_excel(self, path, header, data_list, strt_idx=0):
        workbook = openpyxl.Workbook()
        sheet = workbook.active

        row = 1
        self.make_row(sheet, row, header[strt_idx:])

        for data_arr in data_list:
            row += 1
            self.make_row(sheet, row, data_arr[strt_idx:])

        workbook.save(filename=path + self.ext_xlsx)