#!/usr/bin/env python3

import os
import os.path
import sys
import argparse
import pandas as pd
import numpy as np
import subprocess
from pyfaidx import Fasta, Faidx


def bed_to_fasta(genome, bed, write_path):
    fasta = []
    for line in bed.itertuples():
        chromosome = line[1]
        start = line[2]
        end = line[3]
        ss = genome[chromosome][start:end]
        seq = ss.seq.upper()
        fasta.append(f'>{chromosome}:{start}-{end}\n')
        fasta.append(f'{seq}\n')
    with open(write_path, 'w') as file:
        for line in fasta:
            file.write(line)
    pass


def check_peaks(peaks, idx):
    filter = []
    chromosomes = set(idx.index.keys())
    for line in peaks.itertuples():
        chromosome = line[1]
        start = line[2]
        end = line[3]
        if not chromosome in chromosomes:
            filter.append(False)
            continue
        if end > idx.index[chromosome]['rlen']:
            filter.append(False)
            continue
        filter.append(True)
    filter = np.array(filter)
    filtered_peaks = peaks[np.logical_not(filter)]
    peaks = peaks[filter]
    return peaks, filtered_peaks


def run_dust(fasta_path, thr=20):
    args = ['dust', fasta_path, f'{thr}']
    out = subprocess.run(args, shell=False, capture_output=True)
    fasta = out.stdout.decode()
    fasta = fasta.split('\n')
    container = []
    string = ''
    for i in fasta:
        if i.startswith('>'):
            container.append(string)
            container.append(i)
            string = ''
        else:
            string += i
    container.append(string)
    return container[1:]


def write_fasta(fasta, path):
    with open(path, 'w') as file:
        for line in fasta:
            file.write(f'{line}\n')
    pass


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('output', action='store', help='Output directory to write results')
    parser.add_argument('peaks', action='store', help='Path to PEAKS in BED format')
    parser.add_argument('genome', action='store', help='Path to GENOME in FASTA format')
    parser.add_argument('-c', '--count', action='store', type=int, dest='count',
                        required=False, default=2000, help='The the total number of peaks that will be taken \
                        for the test and train samples from the original data. Only even peaks are taken according to their order for test sample. \
                        Only odd peaks are taken according to their order for train sample.')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()

    write_dir = args.output
    peaks_path = args.peaks
    genome_path = args.genome
    number_of_peaks = args.count

    train_dir = f'{write_dir}/train/'
    test_dir = f'{write_dir}/test/'
    peaks_dir = f'{write_dir}/peaks/'

    if not os.path.exists(write_dir):
        os.mkdir(write_dir)
    if not os.path.exists(train_dir):
        os.mkdir(train_dir)
    if not os.path.exists(test_dir):
        os.mkdir(test_dir)
    if not os.path.exists(peaks_dir):
        os.mkdir(peaks_dir)

    genome = Fasta(genome_path)
    idx = Faidx(genome_path)
    # if not os.path.exists(f'{genome_path}.fai'):
    #     print('FASTA doesn`t have index. Build index for FASTA file')
    #     idx.build(genome_path, f'{genome_path}.fai')


    shift = 100
    peaks = pd.read_csv(peaks_path, sep='\t', header=None, comment='#')
    peaks = peaks.sort_values(by=6, ascending=False)
    peaks.columns = ['chr', 'start',  'end', 'abs_summit', 'pileup', '-log10(pvalue)', 'fold_enrichment', '-log10(qvalue)', 'name', 'supporting_peakcallers']

    new_start = list(peaks['abs_summit'] - shift)
    new_start = [max(i) for i in zip(peaks['start'], new_start)]
    new_end = list(peaks['abs_summit'] + shift)
    new_end = [min(i) for i in zip(peaks['end'], new_end)]

    peaks['start'] = new_start
    peaks['end'] = new_end

    peaks, filtered_peaks = check_peaks(peaks, idx)
    if len(filtered_peaks) > 0:
        number_of_filtered_peaks = len(filtered_peaks)
        print(f'List of peaks ({number_of_filtered_peaks}) filtered due to chromosome length mismatch or incorrect chromosome name.')
        print(filtered_peaks)
    number_of_peaks_actual = len(peaks)

    if number_of_peaks > number_of_peaks_actual:
        print(f'BED file doesn`t have enought peaks.\
                Actual number of peaks is {number_of_peaks_actual}, but {number_of_peaks} is needed')
        if number_of_peaks_actual > 500:
            number_of_peaks = number_of_peaks_actual
            print(f'So only {number_of_peaks_actual} peaks will be used for analysis')
            train_index = np.arange(0, number_of_peaks, 2)
            test_index = np.arange(1, number_of_peaks, 2)
        else:
            print(f'As  number of peaks less than 500. All peaks will be used as train and test')
            train_index = np.arange(0, number_of_peaks_actual)
            test_index = np.arange(0, number_of_peaks_actual)
    else:
        train_index = np.arange(0, number_of_peaks, 2)
        test_index = np.arange(1, number_of_peaks, 2)

    train_bed = peaks.iloc[train_index]
    test_bed = peaks.iloc[test_index]


    train_bed.to_csv(f'{train_dir}/train.bed', header=None, sep='\t', index=False)
    bed_to_fasta(genome, train_bed, f'{train_dir}/train.fa')
    # fasta_masked = run_dust(f'{train_dir}/train.fa', thr=40)
    # write_fasta(fasta_masked, f'{train_dir}/train.fa')

    test_bed.to_csv(f'{test_dir}/test.bed', header=None, sep='\t', index=False)
    bed_to_fasta(genome, test_bed, f'{test_dir}/test.fa')
    # fasta_masked = run_dust(f'{test_dir}/test.fa', thr=40)
    # write_fasta(fasta_masked, f'{test_dir}/test.fa')

    peaks.to_csv(f'{peaks_dir}/peaks.bed', header=None, sep='\t', index=False)
    bed_to_fasta(genome, peaks, f'{peaks_dir}/peaks.fa')
    pass


if __name__ == '__main__':
    main()
