import pickle, sys, os
import argparse
from collections import defaultdict, Counter

import bamnostic as bs


"""following tag names are explained here: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam#bam-align-tags
GN: list of compatible gene name
CB: error-corrected cell barcode
UR: raw UMI
UB: error-corrected UMI
"""

# list of probe index associated with housekeeping genes
genes_to_ignore = [
    '01962',
    '01192',
    '01191',
    '01190',
    '03138',
    '00623',
    '03211',
    '03210',
    '03207',
    '03171',
    '02059',
]

def parse_bam_file(exp_name, data_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    input_file = f'{data_dir}/{exp_name}/possorted_genome_bam.bam'
    output_file = f'{out_dir}/{exp_name}_no_housekeeping.pkl'

    # read bam file
    bam = bs.AlignmentFile(input_file, 'rb')

    # dictionary of counter is indexed by cell barcode and
    # each of this counter is indexed by (umi, probe)-tuple
    #bc_and_raw_umi = defaultdict(Counter)
    bc_and_corrected_umi = defaultdict(Counter)

    n_reads = 0
    n_unmapped_reads = 0
    n_housekeeping_reads = 0

    for read in bam:
        
        n_reads += 1
        
        if gx_field := read.tags.get('GN'):
            probe = gx_field[1]
            probe_num = probe.split('_')[1]
            if probe_num in genes_to_ignore:
                n_housekeeping_reads += 1
                continue
        else:
            n_unmapped_reads += 1
            continue

        if bc_field := read.tags.get('CB'):
            barcode = bc_field[1].split('-')[0]         
            #raw_umi = read.tags.get('UR')[1].split('-')[0]
            #bc_and_raw_umi[barcode][(raw_umi, probe)] += 1
            
            if umi_field := read.tags.get('UB'):
                corrected_umi = umi_field[1].split('-')[0]
                bc_and_corrected_umi[barcode][(corrected_umi, probe)] += 1
       
    result = {
        #'raw_umi': bc_and_raw_umi,
        'corrected_umi': bc_and_corrected_umi,
        'num_total_reads': n_reads,
        'num_housekeeping_reads': n_housekeeping_reads,
        'num_unmapped_reads': n_unmapped_reads,
    }

    pickle.dump(result, open(output_file, 'wb'))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--data_dir', 
        help='path to directory containing experiment folders')
    parser.add_argument(
        '--out_dir', 
        help='path to directory where parsed result will be saved')
    parser.add_argument(
        '--exp_name', 
        help='name of your experiment (=name of the folder in the data_dir')
    args = parser.parse_args()
    parse_bam_file(args.exp_name, args.data_dir, args.out_dir)
