from exalu.run_model import run_ead
from exalu.data_preprocess.curate_data import add_context
import pybedtools
import sys, argparse, textwrap
import os

def get_arguments():
    parser = argparse.ArgumentParser()
    # class SubCommandCommonParser(argparse.ArgumentParser):
    #     def __init__(self, *args, **kwargs):
    #         super(SubCommandCommonParser, self).__init__(*args, **kwargs)
    #         # common options for each sub command
    #         # self.add_
            
    subparsers = parser.add_subparsers(help='select bed or fasta input mode', dest='input_mode')

    parser_bed = subparsers.add_parser('bed', help='infer with bed file')
    parser_bed.add_argument('-b', metavar='ALU_BED_FILE', type=str, required=True, help='the input Alu bed file')
    parser_bed.add_argument('-r', metavar='REF_GENOME_FILE', type=str, required=True, help='the reference genome file, it should be hg38c.fa')
    parser_bed.add_argument('-m', metavar='MODEL_WEIGHTS_FILE', type=str, required=True, help='the trained model weights file')
    parser_bed.add_argument('-o', metavar='OUTPUT_DIR', type=str, required=True, help='the directory contains temp files and final output file')

    parser_fa = subparsers.add_parser('fasta', help='infer with fasta file')
    parser_fa.add_argument('-f', metavar='ALU_FASTA_FILE', type=str, required=True, help='the input Alu fa file')
    parser_fa.add_argument('-m', metavar='MODEL_WEIGHTS_FILE', type=str, required=True, help='the trained model weights file')
    parser_fa.add_argument('-o', metavar='OUTPUT_DIR', type=str, default='./out', help='the directory contains temp files and final output file, default ./out')

    if len(sys.argv) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def infer(model_weights_file, output_dir, alu_bed_file=None, ref_genome_file=None, alu_fa_file=None):
    os.makedirs(output_dir, exist_ok=True)
    if alu_bed_file:
        contexted_bed_file = os.path.join(output_dir, 'temp_contexted_alu.bed')
        fa_file = os.path.join(output_dir, 'temp_contexted_alu.fa')
        ref_genome = pybedtools.example_filename(ref_genome_file)
        add_context(alu_bed_file, contexted_bed_file, 'bothfix', 25)
        alu_bed = pybedtools.BedTool(contexted_bed_file)
        alu_bed.sequence(fi=ref_genome, fo=fa_file, name=True, s=True)
    if alu_fa_file:
        fa_file = alu_fa_file
    run_ead(True, 'cnet', 'bothfix', 25, 4, output_dir, 'simpleinfer', fa_file, 'results.txt', model_weights_file)


if __name__ == "__main__":
    args = get_arguments()
    if args.input_mode == 'bed':
        infer(model_weights_file=args.m, output_dir=args.o, alu_bed_file=args.b, ref_genome_file=args.r)
    if args.input_mode == 'fasta':
        infer(model_weights_file=args.m, output_dir=args.o, alu_fa_file=args.f)
