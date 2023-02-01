from collections import defaultdict
import os
import pybedtools
import csv
import sys
import seaborn
import matplotlib.pyplot as plt


SEQ_LEN_TH_UP = 350
SEQ_LEN_TH_LOW = 100

def split(src_bed_file, len_th):
    '''
    chr22	50622719	50624718	SRR1317771	+	0	chr22	50623919	50624078	h38_mk_AluJ	0	-	159
    '''
    src_bed_fh = open(src_bed_file, 'r')
    src_bed_reader = csv.reader(src_bed_fh, delimiter='\t')
    exon_ids = []
    alu_ids = []
    for idx, row in enumerate(src_bed_reader):
        for i, j in [(1, 4), (7, 11)]:
            l = int(row[i])
            r = int(row[i + 1])
            seq_len = r - l
            if len_th == True and (seq_len > SEQ_LEN_TH_UP or seq_len < SEQ_LEN_TH_LOW):
                continue
            _id = (row[i - 1], l, r, row[j])
            if i == 1:
                exon_ids.append(_id)
            else:
                if _id == ('chr11', '112229324', '112229701', '-'):
                    print(src_bed_file)
                alu_ids.append(_id)
    src_bed_fh.close()
    return exon_ids, alu_ids

def check_dir(dir, len_th):
    unique_set = set()
    dup_num = 0
    # read_file
    file_dict = {}
    sample_counter = 0
    sample_names = set()
    exon_ids, alu_ids = [], []
    for i in os.listdir(dir):
        if i.endswith('_Alu.r.overlap.filtered.bed') or i.endswith('_Alu.r.overlap.filtered.bed.withGeneNames') or \
            i == 'GENCODE.v36.ALUs.all.overlap.filtered.bed' or i == 'rheMac_ENSEMBL_Alu.overlap.filtered.bed':
            sample_counter += 1
            sample_name = i.split('_')[0]
            sample_names.add(sample_name)
            src_bed_file = os.path.join(dir, i)
            sample_exon_ids, sample_alu_ids = split(src_bed_file, len_th)
            exon_ids += sample_exon_ids
            alu_ids += sample_alu_ids
    sample_names_str = ','.join(sample_names)
    print(f'{os.path.basename(dir)}\t{len(sample_names)}\t{sample_names_str}\t{len(alu_ids)}\t{len(set(alu_ids))}', end='\t')
    print(f'{len(exon_ids)}\t{len(set(exon_ids))}')
    return exon_ids, alu_ids

def scan_all_dir(len_th):
    print('Tissue\t#Samples\tSamples\t#Alus\t#Unique_Alus\t#Exons\t#Unique_Exons')
    skip = ['simpleinfer', 'OtherSpecies', 'Upf1SRA', 'MOAT', 'OtherSpecies', 'hnRNPC', 'NMD', 'Gencode', 'Upf1SRA', 'BodyMap']
    dir_lst = []
    work_dir = '/home/zhe29/Projects/eXAlu/data/Fixed_add13'
    exon_ids, alu_ids = [], []
    for dir in os.listdir(work_dir):
        if dir.startswith('padding'):
            continue
        dir_lst.append(dir)
    for dir in dir_lst:
        if dir not in skip:
            dir_exon_ids, dir_alu_ids = check_dir(os.path.join(work_dir, dir), len_th)
            exon_ids += dir_exon_ids
            alu_ids += dir_alu_ids
    print(f'ALL\tN/A\tN/A\t{len(alu_ids)}\t{len(set(alu_ids))}', end='\t')
    print(f'{len(exon_ids)}\t{len(set(exon_ids))}')

    if len_th == True:
        return
    # draw the distribution
    unique_exon_ids = list(set(exon_ids))
    unique_alu_ids = list(set(alu_ids))
    len_dict = defaultdict(list)

    for k, ids in [('alu', alu_ids), ('unique_alu', unique_alu_ids)]:
    # for k, ids in [('exon', exon_ids), ('unique_exon', unique_exon_ids), ('alu', alu_ids), ('unique_alu', unique_alu_ids)]:
        for _id in ids:
            l = _id[2] - _id[1]
            len_dict[k].append(l)
        seaborn.histplot(data=len_dict[k], binwidth=10)
        plt.savefig(f'./data/output/imgs/data_len_distribution_{k}_abstract.png')
        plt.close()


if __name__ == '__main__':
    scan_all_dir(len_th=True)
