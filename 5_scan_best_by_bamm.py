#!/usr/bin/env python3

import re
import subprocess
import os
import sys
import shutil
import random
import bisect
import argparse
import operator
import itertools
from math import log
from shutil import copyfile
from operator import itemgetter


# BAMM MODEL
def parse_bamm_and_bg_from_file(bamm_file, bg_file):

    # Read BaMM file
    if os.path.isfile(bamm_file):
        motif_order = 0
        with open(bamm_file) as file:
            for line in file:
                if line[0] != '\n':
                    motif_order = motif_order + 1
                else:
                    break

        # count the motif length
        motif_length = int(sum(1 for line in open(bamm_file)) / (motif_order + 1))

        # read in bamm model
        model = {}
        for k in range(motif_order):
            model[k] = []

        with open(bamm_file) as file:
            for j in range(motif_length):
                for k in range(motif_order):
                    model[k].append([float(p) for p in file.readline().split()])
                file.readline()

    else:
        print('File {0} does not exist'.format(bg_file))
        sys.exit()

    # Read BG file
    bg = {}
    order = 0
    if os.path.isfile(bg_file):
        with open(bg_file) as bgmodel_file:
            line = bgmodel_file.readline()  # skip the first line for K
            line = bgmodel_file.readline()  # skip the second line for Alpha
            while order < motif_order:
                line = bgmodel_file.readline()
                bg_freq = [float(p) for p in line.split()]
                bg[order] = bg_freq
                order += 1
    else:
        print('File {0} does not exist'.format(bg_file))
        sys.exit()
    return(model, bg, order-1)


def make_k_mers(order):
    #  make list with possible k-mer based on bHMM model
    tmp = itertools.product('ACGT', repeat=order + 1)
    k_mer = []
    for i in tmp:
        k_mer.append(''.join(i))
    k_mer_dict = dict()
    index = 0
    for i in k_mer:
        k_mer_dict[i] = index
        index += 1
    return(k_mer_dict)


def make_log_odds_bamm(bamm, bg):
    log_odds_bamm = dict()
    for order in bamm.keys():
        log_odds_bamm[order] = [list(map(lambda x: log(x[0] / x[1], 2), zip(bamm_col, bg[order]))) for bamm_col in bamm[order]]
    return(log_odds_bamm)


def bamm_to_dict(log_odds_bamm, order, k_mers):
    bamm_dict = {}
    for k_mer in k_mers:
        bamm_dict[k_mer] = list()
    for i in range(len(log_odds_bamm[order])):
        for index, k_mer in enumerate(k_mers):
            bamm_dict[k_mer].append(log_odds_bamm[order][i][index])
    return(bamm_dict)


def read_bamm(bamm_path, bg_path):
    bamm, bg, order = parse_bamm_and_bg_from_file(bamm_path, bg_path)
    container = dict()
    log_odds_bamm = make_log_odds_bamm(bamm, bg)
    for i in range(order + 1):
        k_mers = make_k_mers(i)
        bamm_dict = bamm_to_dict(log_odds_bamm, i, k_mers)
        container.update(bamm_dict)
    return(container, order)


def score_bamm(site, bamm, order, length_of_site):
    score = 0
    for index in range(order):
        score += bamm[site[0:index + 1]][index]
    for index in range(length_of_site - order):
        score += bamm[site[index:index+order + 1]][index + order]
    return(score)


def min_score_bamm(bamm, order, length_of_site):
    scores = []
    min_score = 0
    for index in range(order):
        k_mers = itertools.product('ACGT', repeat=index + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for index in range(length_of_site - order):
        k_mers = itertools.product('ACGT', repeat=order + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for s in scores:
        min_score += min(s)
    return(min_score)


def max_score_bamm(bamm, order, length_of_site):
    scores = []
    max_score = 0
    for index in range(order):
        k_mers = itertools.product('ACGT', repeat=index + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for index in range(length_of_site - order):
        k_mers = itertools.product('ACGT', repeat=order + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for s in scores:
        max_score += max(s)
    return(max_score)


def to_score(norm_value, bamm, order, length_of_site):
    max_s = max_score_bamm(bamm, order, length_of_site)
    min_s = min_score_bamm(bamm, order, length_of_site)
    score = norm_value * (max_s - min_s) + min_s
    return(score)


def to_norm(score, bamm, order, length_of_site):
    max_s = max_score_bamm(bamm, order, length_of_site)
    min_s = min_score_bamm(bamm, order, length_of_site)
    norm_value = (score - min_s) / (max_s - min_s)
    return(norm_value)


def read_fasta(path):
    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            #print(line)
            if line.startswith('>'):
                line = line[1:].strip().split()
                record = dict()
                record['name'] = line[0]
            else:
                record['seq'] = line.strip().upper()
                fasta.append(record)
    file.close()
    return(fasta)


def complement(record):
    output = dict(record)
    seq = output['seq'].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]
    output['seq'] = seq
    return(output)


def check_nucleotides(site):
    s = set(site)
    n = {'A', 'C', 'G', 'T'}
    if len(s - n) == 0:
        return(True)
    else:
        return(False)


def scan_seqs_by_bamm(record, log_odds_bamm, order):
    motif_length = len(log_odds_bamm[list(log_odds_bamm.keys())[0]])
    reverse_record = complement(record)
    seq = record['seq']
    reverse_seq = reverse_record['seq']
    results = []
    threshold = -1000000
    # scan first strand
    for i in range(len(seq) - motif_length + 1):
        site_seq = seq[i:motif_length + i]
        if not check_nucleotides(site_seq):
            continue
        s = score_bamm(site_seq, log_odds_bamm, order, motif_length)
        if s >= threshold:
            site_dict = dict()
            site_dict['name'] = record['name']
            site_dict['score'] = s
            threshold = s
    # scan second strand
    for i in range(len(seq) - motif_length + 1):
        site_seq = reverse_seq[i:motif_length + i]
        if not check_nucleotides(site_seq):
            continue
        s = score_bamm(site_seq, log_odds_bamm, order, motif_length)
        if s >= threshold:
            site_dict = dict()
            site_dict['name'] = record['name']
            site_dict['score'] = s
            threshold = s
    results.append(site_dict)
    return(results)


def write_list(path, data):
    scores = [i['score'] for i in data]
    with open(path, "w") as file:
        for line in scores:
            file.write("{0}\n".format(line))
    file.close()
    return(0)


def write_results(path, data, tf, log_odds_bamm, order, motif_length):
    with open(path, "w") as file:
        file.write(f'tag\t{tf}\n')
        for line in data:
            tag = line['name']
            score = to_norm(line['score'], log_odds_bamm, order, motif_length)
            file.write(f'{tag}\t{score:.5f}\n')
    pass


def scan_best_by_bamm(results_path, bamm_dir, fasta_path, tf):

    bamm_path = f'{bamm_dir}/bamm_model.ihbcp'
    bg_path = f'{bamm_dir}/bamm.hbcp'

    fasta = read_fasta(fasta_path)
    log_odds_bamm, order = read_bamm(bamm_path, bg_path)
    motif_length = len(log_odds_bamm[list(log_odds_bamm.keys())[0]])
    results = [scan_seqs_by_bamm(record=record, log_odds_bamm=log_odds_bamm, order=order) for record in fasta]
    results = [i for i in results if i != []]
    results = [j for sub in results for j in sub]
    write_results(results_path, results, tf, log_odds_bamm, order, motif_length)
    return(0)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('bamm', action='store', help='Dir with BaMM model')
    parser.add_argument('fasta', action='store', help='path to FSATA')
    parser.add_argument('output', action='store', help='Path to write results')
    parser.add_argument('tf', action='store', help='Name of tf')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()

    results_path = args.output
    bamm_dir = args.bamm
    fasta_path = args.fasta
    tf = args.tf

    scan_best_by_bamm(results_path, bamm_dir, fasta_path, tf)

if __name__ == '__main__':
    main()
