import os
import pybedtools
from .read_tissue_alu import add_context, read_tissue_alu_context

# data_dir = os.getcwd() + '/data'
# fixed_dir = data_dir + '/Fixed_add13/' # stable model use this one
# fixed_dir = data_dir + '/Fixed_genesplicer/' # just used for generating the gene splicer analysis

# human_alu_bed_file = data_dir + '/Annotation/Alu_repeats_hg38_filtered.bed'

def curate_each_tissue(genome, strand, work_dir, tissue):
    '''
    This function reads data from each tissue's directory (src_pos_bed_tissue_dir), and combine all samples within a same directory below
    f'padding_{mode}_{single_side_pad_len}_{testsets}/{tissue}'
    '''
    # read_file
    src_pos_bed_tissue_dir = os.path.join(work_dir, 'after_read_Fixed', tissue, 'split_bed')

    dir_str = f'{work_dir}/{tissue}'
    dst_pos_alu_fa = dir_str + '_pos_alu.fa' # work_dir/tissue_pos_alu.fa
    dst_pos_exon_fa = dir_str + '_pos_exon.fa'
    dst_pos_alu_bed = dir_str + '_pos_alu.bed'
    dst_pos_exon_bed = dir_str + '_pos_exon.bed'

    ref_fa = pybedtools.example_filename(genome)
    alu_bed = None
    exon_bed = None
    # get input files prefix list
    in_file_prefix_lst = []
    for i in os.listdir(src_pos_bed_tissue_dir):
        prefix = i.split('_')[0]
        if prefix not in in_file_prefix_lst:
            in_file_prefix_lst.append(prefix)
    dst_dir_name = os.path.dirname(dst_pos_alu_fa)
    if not os.path.isdir(dst_dir_name):
        os.makedirs(dst_dir_name)
    with open(dst_pos_alu_bed,  'w') as out_alu_file, \
         open(dst_pos_exon_bed, 'w') as out_exon_file:
        for i in os.listdir(src_pos_bed_tissue_dir):
            # zero_file_prefix = i.split('_')[0]
            if i.endswith('alu.bed'):
                with open(os.path.join(src_pos_bed_tissue_dir, i), 'r') as in_alu_file:
                    for line in in_alu_file:
                        out_alu_file.write(line)
            if i.endswith('exon.bed'):
                with open(os.path.join(src_pos_bed_tissue_dir, i), 'r') as in_exon_file:
                    for line in in_exon_file:
                        out_exon_file.write(line)
    alu_bed = pybedtools.BedTool(dst_pos_alu_bed)
    exon_bed = pybedtools.BedTool(dst_pos_exon_bed)
    alu_bed.sequence(fi=ref_fa, fo=dst_pos_alu_fa, name=True, s=strand)
    exon_bed.sequence(fi=ref_fa, fo=dst_pos_exon_fa, name=True, s=strand)


def curate_whole(genome, strand, work_dir, mode, single_side_pad_len, tissue_lst, human_alu_bed_file, infer_set_lst=['Gencode', 'MOAT']):
    '''
    This function curates whole dataset from every tissues
    Gencode is used for inference, so it should not be included in trainset
    pos: all exclude Gencode
    neg: all_alu - all_pos_alu, this all_pos_alu include Gencode,
        because we don't want infer set and train set share same neg data.
    
    data/curated_data/padding____/whole_pos_alu.fa doesn't have __infer sets__
    data/curated_data/padding____/whole_pos_alu.bed has __infer sets__

    '''
    whole_dir_str = f'{work_dir}/whole'
    dst_whole_pos_alu_bed  = whole_dir_str + '_pos_alu.bed'
    dst_whole_pos_exon_bed = whole_dir_str + '_pos_exon.bed'
    dst_whole_pos_alu_fa   = whole_dir_str + '_pos_alu.fa'
    dst_whole_pos_exon_fa  = whole_dir_str + '_pos_exon.fa'
    alu_bed = None
    exon_bed = None
    ref_fa = pybedtools.example_filename(genome)
    out_alu_file = open(dst_whole_pos_alu_bed, 'w')
    out_exon_file = open(dst_whole_pos_exon_bed, 'w')
    for tissue in tissue_lst:
        if tissue in infer_set_lst:
            # skip Gencode and MOAT; developing;
            continue # skip
        tissue_dir_str = f'{work_dir}/{tissue}'
        src_tissue_pos_alu_bed  = tissue_dir_str + '_pos_alu.bed'
        src_tissue_pos_exon_bed = tissue_dir_str + '_pos_exon.bed'
        print('adding', tissue_dir_str)
        with open(src_tissue_pos_alu_bed, 'r') as in_alu_file:
            for line in in_alu_file:
                out_alu_file.write(line)
        with open(src_tissue_pos_exon_bed, 'r') as in_exon_file:
            for line in in_exon_file:
                out_exon_file.write(line)
    out_alu_file.close()
    out_exon_file.close()

    alu_bed = pybedtools.BedTool(dst_whole_pos_alu_bed)
    exon_bed = pybedtools.BedTool(dst_whole_pos_alu_bed)
    # this is the fasta for training, not including infer set
    alu_bed.sequence(fi=ref_fa, fo=dst_whole_pos_alu_fa, name=True, s=strand)
    exon_bed.sequence(fi=ref_fa, fo=dst_whole_pos_exon_fa, name=True, s=strand)

    # add infer set for preparing neg data
    for infer_name in infer_set_lst:
        if infer_name == 'OtherSpecies':
            break
        if infer_name in tissue_lst:
            out_alu_file = open(dst_whole_pos_alu_bed, 'a')
            out_exon_file = open(dst_whole_pos_exon_bed, 'a')
            infer_dir_str = f'{work_dir}/{infer_name}'
            print('adding', infer_dir_str)
            src_infer_pos_alu_bed  = infer_dir_str + '_pos_alu.bed'
            src_infer_pos_exon_bed = infer_dir_str + '_pos_exon.bed'
            with open(src_infer_pos_alu_bed, 'r') as in_alu_file:
                for line in in_alu_file:
                    out_alu_file.write(line)
            with open(src_infer_pos_exon_bed, 'r') as in_exon_file:
                for line in in_exon_file:
                    out_exon_file.write(line)
            out_alu_file.close()
            out_exon_file.close() 
    alu_bed = pybedtools.BedTool(dst_whole_pos_alu_bed)
    exon_bed = pybedtools.BedTool(dst_whole_pos_alu_bed)

    # methond 1. merge alu and exon to get all bed file, TODO:
    # then use human alu to substract the all bed file.
    # with open(pos_all_bed_file.format(whole_dir), 'w') as out_file:
        # for f in [pos_alu_bed_file.format(whole_dir), pos_exon_bed_file.format(whole_dir)]:
            # with open(f, 'r') as in_file:
                # for line in in_file:
                    # out_file.write(line)
    # pos_all_bed = pybedtools.BedTool(pos_all_bed_file.format(whole_dir))
    # pos_all_bed.sequence(fi=ref_fa, fo=pos_all_fa_file.format(whole_dir), name=True, s=strand)
    # human_alu_bed = pybedtools.BedTool(human_alu_bed_file)
    # neg_alu_bed = human_alu_bed.subtract(pos_all_bed, A=True)

    # method 2. use human alu to substract the alu bed file, without exon file
    human_alu_bed_file_context = os.path.join(work_dir, os.path.basename(human_alu_bed_file))
    add_context(human_alu_bed_file, human_alu_bed_file_context, mode=mode, single_side_pad_len=single_side_pad_len)
    human_alu_bed = pybedtools.BedTool(human_alu_bed_file_context)
    print('Neg Alu candidates for training:', human_alu_bed_file)
    neg_alu_bed = human_alu_bed.subtract(alu_bed, A=True)
    print('Length before substract pos alus: ', human_alu_bed.count())
    print('Length before substract pos alus: ', neg_alu_bed.count())
    # get neg sequence
    dst_neg_bed = os.path.join(work_dir, 'neg_alu.bed')
    dst_neg_fa  = os.path.join(work_dir, 'neg_alu.fa')
    neg_alu_bed.saveas(dst_neg_bed) # bed file, not needed
    neg_alu_bed.sequence(fi=ref_fa, fo=dst_neg_fa, name=True, s=strand)

def curate_whole_otherspecies(genome, strand, work_dir, mode, single_side_pad_len, otherspecies='OtherSpecies'):
    '''
    deprecated;
    This function curates whole dataset from every tissues
    Gencode is used for inference, so it should not be included in trainset
    pos: all exclude Gencode
    neg: all_alu - all_pos_alu, this all_pos_alu include Gencode,
        because we don't want infer set and train set share same neg data.
    
    data/curated_data/padding____/whole_pos_alu.fa doesn't have __infer sets__
    data/curated_data/padding____/whole_pos_alu.bed has __infer sets__

    '''
    alu_bed = None
    exon_bed = None
    ref_fa = pybedtools.example_filename(genome)
    tissue_dir_str = f'{work_dir}/{otherspecies}'
    src_tissue_pos_alu_bed  = tissue_dir_str + '_pos_alu.bed'
    src_tissue_pos_exon_bed = tissue_dir_str + '_pos_exon.bed'
    dst_whole_pos_alu_fa   = tissue_dir_str + '_pos_alu.fa'
    dst_whole_pos_exon_fa  = tissue_dir_str + '_pos_exon.fa'
    print('adding', tissue_dir_str)
    alu_bed = pybedtools.BedTool(src_tissue_pos_alu_bed)
    exon_bed = pybedtools.BedTool(src_tissue_pos_exon_bed)
    alu_bed.sequence(fi=ref_fa, fo=dst_whole_pos_alu_fa, name=True, s=strand)
    exon_bed.sequence(fi=ref_fa, fo=dst_whole_pos_exon_fa, name=True, s=strand)

    # method 2. use human alu to substract the alu bed file, without exon file
    # TODO: hl
    otherspecies_alu_bed = pybedtools.BedTool('/home/zitong/Projects/eXAlu/data/OtherSpecies/rheMac10_Alu_annotations.bed')
    neg_alu_bed = otherspecies_alu_bed.subtract(alu_bed, A=True)
    # get neg sequence
    neg_dir_str = f'{work_dir}/{otherspecies}'
    # TODO:
    # can shorten below
    # TODO:
    # neg alu wo context now
    # check curate_whole()
    dst_neg_bed = neg_dir_str + f'_neg_alu.bed'
    dst_neg_fa  = neg_dir_str + f'_neg_alu.fa'
    neg_alu_bed.saveas(dst_neg_bed) # bed file, not needed
    neg_alu_bed.sequence(fi=ref_fa, fo=dst_neg_fa, name=True, s=strand)

def curate_simpleinfer(strand, mode, single_side_pad_len, infer_bed_file=None, work_dir=None, data_dir=None):
    genome = data_dir + '/shared/hg38/hg38c.fa' 
    # read_tissue_alu_context(genome=genome, mode=mode, single_side_pad_len=single_side_pad_len, tissue='simpleinfer')
    ref_fa = pybedtools.example_filename(genome)
    bed_name_wo_ext = os.path.basename(infer_bed_file).split('.')[0]
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    dst_bed = os.path.join(work_dir, bed_name_wo_ext + '.bed')
    dst_fa = os.path.join(work_dir, bed_name_wo_ext + '.fa')
    if infer_bed_file != None:
        add_context(src_bed_file=infer_bed_file, dst_bed_file=dst_bed, mode=mode, single_side_pad_len=single_side_pad_len)
    print('adding', dst_fa)
    alu_bed = pybedtools.BedTool(dst_bed)
    alu_bed.sequence(fi=ref_fa, fo=dst_fa, name=True, s=strand)

def curate_data(infer_set, strand, mode='none', single_side_pad_len=0, work_dir=None, data_dir=None):
    '''
    mode: left - only left context; right - only right context; both - left and right context
        none - no context
    '''
    print(work_dir)
    tissue_lst = []
    # dataset (tissue or source) selection
    # fixed_dir = data_dir + '/Fixed_add13/' # stable model use this one
    training_bed_dir = data_dir +  '/Fixed_add13/'
    for d in os.listdir(training_bed_dir):
        if d.startswith('padding_') or d in ['f0', 'f350', 'simpleinfer']:
            continue
        if infer_set == 'MOAT':
            if d in ['OtherSpecies', 'hnRNPC', 'NMD', 'Gencode', 'Upf1SRA', 'BodyMap']:
                print('skipping ..............', d)
                continue
        if infer_set == 'Gencode':
            if d in ['MOAT', 'OtherSpecies', 'hrRNPC', 'NMD']:
                continue
        if infer_set == 'OtherSpecies':
            if d in ['Gencode', 'MOAT']:
                continue
        if infer_set == 'all_human_analysis':
            if d in ['OtherSpecies']:
                continue
        tissue_lst.append(d)
    print(tissue_lst)
    for tissue in tissue_lst:
        genome = data_dir + '/shared/hg38/hg38c.fa'
        read_tissue_alu_context(genome, strand, mode, single_side_pad_len, tissue, training_bed_dir, work_dir)
        curate_each_tissue(genome, strand, work_dir, tissue)
    genome = data_dir + '/shared/hg38/hg38c.fa' 

    # old version 04/09/22 Alu_repeats_hg38_filtered.bed may still have false negative data
    human_alu_bed_file = data_dir + '/Annotation/Alu_repeats_hg38_filtered.bed'
    # new version 10/21/22 intergenic alus: Alu_repeats_hg38_filtered.bed - gencode_exons - gencode_introns
    # human_alu_bed_file = data_dir + '/Annotation/intergenic_alu/all_alu_minus_exon_minus_intron.bed'
    curate_whole(genome, strand, work_dir, mode, single_side_pad_len, tissue_lst, human_alu_bed_file=human_alu_bed_file,infer_set_lst=[infer_set])
