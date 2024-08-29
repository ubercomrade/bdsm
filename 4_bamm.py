#!/usr/bin/env python3

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


def calculate_prc(true_score, false_score):
    table = []
    prc = {'PRECISION': [0], 'RECALL': [0]}
    number_of_yes = len(true_score)
    yes = 0
    no = 0
    for t in true_score:
        table.append({'score': t, 'type': 1})
    for f in false_score:
        table.append({'score': f, 'type': 0})
    table.sort(key=operator.itemgetter('score'), reverse=True)
    score = table[0]['score']
    for line in table:
        if score != line['score']:
            score = line['score']
            precision = yes / (yes + no)
            recall = yes / number_of_yes
            prc['PRECISION'].append(precision)
            prc['RECALL'].append(recall)
        if line['type'] == 1:
            yes +=1
        else:
            no +=1
    precision = yes / (yes + no)
    recall = yes / number_of_yes
    prc['PRECISION'].append(precision)
    prc['RECALL'].append(recall)
    prc['PRECISION'][0] = prc['PRECISION'][1]
    return prc


def calculate_fprs(true_scores, false_scores):
    fprs = []
    false_scores.sort()
    number_of_sites = len(false_scores)
    #true_scores_uniq = list(set(true_scores))
    #true_scores_uniq.sort(reverse=True)
    true_scores.sort(reverse=True)
    for score in true_scores:
        fpr = (number_of_sites - bisect.bisect_right(false_scores, score)) / number_of_sites
        if fpr == 0:
            fprs.append(0.5 / number_of_sites)
        else:
            fprs.append(fpr)
    return(fprs)


def write_roc(path, data):
    header = list(data.keys())
    with open(path, 'w') as file:
        file.write('\t'.join(header) + '\n')
        for i in zip(*data.values()):
            file.write('\t'.join((map(str,i))) + '\n')
    return(0)


def calculate_short_roc(fprs, step=1):
    table = {'TPR': [0], 'FPR': [0]}#, 'SITES': [0]}
    current_number_of_sites = 0
    total_number_of_sites = len(fprs)
    for i in range(1, 100, step):
        position = round(total_number_of_sites * (i / 100))
        table['TPR'].append(i / 100)
        table['FPR'].append(fprs[position - 1])
        #table['SITES'].append(len(fprs[:position]))
    return(table)


def calculate_staged_roc(fprs, step=1):
    table = {'TPR': [0], 'FPR': [0]}#, 'SITES': [0]}
    current_number_of_sites = 0
    total_number_of_sites = len(fprs)
    for i in range(1, 100, step):
        position = round(total_number_of_sites * (i / 100))
        table['TPR'].append(i / 100)
        table['FPR'].append(fprs[position - 1])
        #table['SITES'].append(len(fprs[:position]))
    return(table)


def write_auc(path, auc, length):
    with open(path, 'a') as file:
        file.write('{0}\t{1}\n'.format(length, auc))
    pass


def write_auc_with_order(path, auc, length, order):
    with open(path, 'a') as file:
        file.write('{0}\t{1}\t{2}\n'.format(order, length, auc))
    pass

def write_table(path, *args):
    with open(path, 'a') as file:
        file.write('\t'.join(list(map(str, args))[::-1])+'\n')
    return(0)

def calculate_merged_roc(fprs):
    table = {'TPR': [0], 'FPR': [0], 'SITES': [0]}
    current_number_of_sites = 0
    total_number_of_sites = len(fprs)
    for i in range(0, total_number_of_sites):
        if fprs[i] > fprs[i - 1]:
            table['TPR'].append(current_number_of_sites / total_number_of_sites)
            table['FPR'].append(fprs[i - 1])
            table['SITES'].append(current_number_of_sites)
        current_number_of_sites += 1
    table['TPR'].append(current_number_of_sites / total_number_of_sites)
    table['FPR'].append(fprs[i])
    table['SITES'].append(current_number_of_sites)
    return(table)

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


def make_pcm(motifs):
    matrix = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        matrix[i[0]] = []
    len_of_motif = len(motifs[0])
    for i in matrix.keys():
        matrix[i] = [0]*len_of_motif
    for i in range(len_of_motif):
        for l in motifs:
            matrix[l[i]][i] += 1
    return(matrix)


def make_pfm(pcm):
    number_of_sites = [0] * len(pcm['A'])
    for key in pcm.keys():
        for i in range(len(pcm[key])):
            number_of_sites[i] += pcm[key][i]
    pfm = dict()
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pfm[i[0]] = []
    first_key = list(pcm.keys())[0]
    nuc_pseudo = 1/len(pcm.keys())
    for i in range(len(pcm[first_key])):
        for nuc in pcm.keys():
            pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites[i] + 1))
    return(pfm)


def read_peaks(path):
    container = []
    letters = {'A', 'C', 'G', 'T'}
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                line = ''.join([l if l in letters else 'N' for l in line.strip().upper()])
                append(line.strip().upper())
    return(container)


def write_fasta(peaks, path):
    with open(path, 'w') as file:
        for index, p in enumerate(peaks):
            file.write('>{}\n'.format(index))
            file.write(p + '\n')
    return(0)


def read_fasta(path):
    container = []
    letters = {'A', 'C', 'G', 'T'}
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                line = ''.join([l if l in letters else 'N' for l in line.strip().upper()])
                append(line.strip().upper())
    return(container)


def write_meme(output, tag, pfm, background, nsites):
    with open(output + '/' + tag + '.meme', 'w') as file:
        file.write('MEME version 4\n\nALPHABET= ACGT\n\nBackground letter frequencies\n')
        file.write('A {0} C {1} G {2} T {3}\n\n'.format(background['A'], background['C'],
                                                        background['G'], background['T']))
        file.write('MOTIF {0}\n'.format(tag))
        file.write(
            'letter-probability matrix: alength= 4 w= {0} nsites= {1}\n'.format(len(pfm['A']), nsites))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0:.8f}\t{1:.8f}\t{2:.8f}\t{3:.8f}\n'.format(i[0], i[1], i[2], i[3]))


def read_seqs_with_complement(path):
    container = []
    letters = {'A', 'C', 'G', 'T'}
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                line = ''.join([l if l in letters else 'N' for l in line.strip().upper()])
                append(line.strip().upper())
                append(complement(line.strip().upper()))
    return(container)


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def creat_background(peaks, length_of_site, counter):
    shuffled_peaks = []
    number_of_sites = 0
    while counter > number_of_sites:
        peak = random.choice(peaks)
        shuffled_peak = ''.join(random.sample(peak, len(peak)))
        shuffled_peaks.append(shuffled_peak)
        number_of_sites += (len(''.join(shuffled_peak)) - length_of_site + 1) * 2
    return(shuffled_peaks)


def calculate_particial_auc(tprs, fprs, pfpr):
    auc = 0
    tpr_old = tprs[0]
    fpr_old = fprs[0]
    for tpr_new, fpr_new in zip(tprs[1:], fprs[1:]):
        if fpr_new >= pfpr:
            auc += (tpr_new + tpr_old) * ((pfpr - fpr_old) / 2)
            break
        auc += (tpr_new + tpr_old) * ((fpr_new - fpr_old) / 2)
        fpr_old = fpr_new
        tpr_old = tpr_new
    if fpr_new < pfpr:
        auc += (2* tpr_new) * ((pfpr - fpr_new) / 2)
    return(auc)


def run_streme(fasta_path, backgroud_path, dir_out, motif_length):
    args = ['streme', '--p', fasta_path,
            '--n', backgroud_path,
           '--oc', dir_out,
           '--objfun', 'de',
           '--w', str(motif_length),
           '-nmotifs', '5']
    p = subprocess.run(args, shell=False, capture_output=True)
    return(0)


def run_streme_hmm_background(fasta_path, dir_out, motif_length):
    args = ['streme', '--p', fasta_path,
            '--kmer', '4',
           '--oc', dir_out,
           '--objfun', 'de',
           '--w', str(motif_length),
           '-nmotifs', '5']
    p = subprocess.run(args, shell=False, capture_output=True)
    return(0)


def parse_streme(path):
    pfm = {'A': [], 'C': [], 'G': [], 'T': []}
    with open(path) as file:
        for line in file:
            if line.startswith('Background'):
                line = file.readline().strip().split()
                background = {'A': float(line[1]),
                             'C': float(line[3]),
                             'G': float(line[5]),
                             'T': float(line[7])}
            elif line.startswith('letter-probability'):
                break
        line = line.strip().split()
        length = int(line[5])
        nsites = int(line[7])
        for i in range(length):
            line = file.readline().strip().split()
            for letter, value in zip(pfm.keys(), line):
                    pfm[letter].append(float(value))
    file.close()
    return pfm, background, length, nsites


def create_bamm_model(peaks_path, backgroud_path, directory, order, meme, extend, basename):
    args = ['BaMMmotif', directory, peaks_path, '--PWMFile', meme,
            '--EM', '--order', str(order), '--Order', str(order),
           '--extend', str(extend), '--basename', str(basename),
           '--negSeqFile', backgroud_path]
    #print(' '.join(args))
    r = subprocess.run(args, capture_output=True)
    bamm_path = directory + '/{}_motif_1.ihbcp'.format(basename)
    bg_path = directory + '/{}.hbcp'.format(basename)
    bamm, order = read_bamm(bamm_path, bg_path)
    return(bamm, order)



def false_scores_bamm(peaks, bamm, order, length_of_site):
    false_scores = []
    append = false_scores.append
    for peak in peaks:
        complement_peak = complement(peak)
        full_peak = peak + 'N' * length_of_site + complement_peak
        n = len(full_peak) - length_of_site + 1
        for i in range(n):
            site = full_peak[i:length_of_site + i]
            if not len(set(site) - {'A', 'C', 'G', 'T'}) == 0:
                continue
            score = score_bamm(site, bamm, order, length_of_site)
            false_scores.append(score)
    return(false_scores)


def best_scores_bamm(peaks, bamm, order, length_of_site):
    true_scores = []
    for peak in peaks:
        complement_peak = complement(peak)
        best = -1000000
        full_peak = peak + 'N' * length_of_site + complement_peak
        n = len(full_peak) - length_of_site + 1
        for i in range(n):
            site = full_peak[i:length_of_site + i]
            if not len(set(site) - {'A', 'C', 'G', 'T'}) == 0:
                continue
            score = score_bamm(site, bamm, order, length_of_site)
            if score >= best:
                best = score
        true_scores.append(best)
    return(true_scores)


def fpr_at_tpr(true_scores, false_scores, tpr):
    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    score = true_scores[round(true_length * tpr) - 1]
    actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
    fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
    return(fpr)


def learn_optimized_bamm_support(train_path, test_path, backgroud_path,
            tmp_dir, write_dir, pwm_dir,counter, order, length, fpr_threshold):

    test = read_fasta(test_path)
    meme = pwm_dir + f'/motif_{length}.meme'

    if os.path.isfile(backgroud_path):
        backgroud = read_peaks(backgroud_path)
    else:
        backgroud = creat_background(tes, length, counter)
        background_path = tmp_dir + '/background.fasta'
        write_fasta(shuffled, background_path)

    if not os.path.isfile(meme):
        run_streme(train_path, backgroud_path, pwm_dir, length)

    # de novo bamm
    bamm, order = create_bamm_model(train_path, backgroud_path, tmp_dir, order, meme, 0, length)

    shutil.copy(tmp_dir + '/{}_motif_1.ihbcp'.format(length),
           f'{write_dir}/statistics/' + f'/bamm_model_{order}_{length}.ihbcp')
    shutil.copy(tmp_dir + '/{}.hbcp'.format(length),
           f'{write_dir}/statistics/' + f'/bamm_{order}_{length}.hbcp')

    true_scores = best_scores_bamm(test, bamm, order, length)
    false_scores = best_scores_bamm(backgroud, bamm, order, length)

    #Calculate FPRs -> ROC -> AUC
    fprs = calculate_fprs(true_scores, false_scores)
    merged_roc = calculate_merged_roc(fprs)
    auc = calculate_particial_auc(merged_roc['TPR'], merged_roc['FPR'], 1.0)

    prc = calculate_prc(true_scores, false_scores)
    auc_prc = calculate_particial_auc(prc['PRECISION'], prc['RECALL'], 1.0)

    # vs frequencies of sites
    false_scores = false_scores_bamm(backgroud, bamm, order, length)
    fprs = calculate_fprs(true_scores, false_scores)
    staged_roc = calculate_staged_roc(fprs, step=1)
    merged_roc = calculate_merged_roc(fprs)
    pauc = calculate_particial_auc(merged_roc['TPR'], merged_roc['FPR'], fpr_threshold)

    #Write ROC and pAUC
    print(f"Length {length}; Order {order}; pAUC at {fpr_threshold} = {pauc}; aucROC = {auc:.4f}; aucPRC = {auc_prc:.4f}")
    write_table(f'{write_dir}/statistics/auc.txt', (pauc / 0.001 + auc + auc_prc), length, order)
    write_roc(f'{write_dir}/statistics/roc_merged_{length}.txt', merged_roc)
    write_roc(f'{write_dir}/statistics/roc_{length}.txt', staged_roc)
    shutil.rmtree(tmp_dir)
    return(0)


def learn_optimized_bamm(train_path, test_path, backgroud_path, tmp_dir, write_dir, pwm_dir,
                        fpr_threshold=0.001, min_length=8, max_length=20, step_length=4, min_order=1, max_order=5):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    else:
        shutil.rmtree(tmp_dir)
        os.mkdir(tmp_dir)
    if not os.path.isdir(write_dir):
        os.mkdir(write_dir)
    if not os.path.isdir(f'{write_dir}/statistics'):
        os.mkdir(f'{write_dir}/statistics')
    if os.path.exists(f'{write_dir}/statistics/auc.txt'):
        os.remove(f'{write_dir}/statistics/auc.txt')
    if not os.path.isdir(pwm_dir):
        os.mkdir(pwm_dir)
    counter = 1000000
    # if os.path.exists(output_auc + '/statistics.txt'):
    #     os.remove(output_auc + '/statistics.txt')
    for order in range(min_order, max_order + 1):
        for length in range(min_length, max_length + 1, step_length):
            run = learn_optimized_bamm_support(train_path, test_path, backgroud_path, tmp_dir, write_dir, pwm_dir,counter, order, length, fpr_threshold)


    length, order = choose_best_model(f'{write_dir}/statistics/auc.txt')
    print(f'Model have best perfomance with length = {length} and order = {order}')
    shutil.copy(f'{pwm_dir}/motif_{length}.meme',
           f'{write_dir}/pwm_model.meme')
    shutil.copy(f'{write_dir}/statistics/bamm_model_{order}_{length}.ihbcp',
           f'{write_dir}/bamm_model.ihbcp')
    shutil.copy(f'{write_dir}/statistics/bamm_{order}_{length}.hbcp',
           f'{write_dir}/bamm.hbcp')
    pass


def choose_best_model(output_auc):
    auc = []
    with open(output_auc, 'r') as file:
        for line in file:
            auc.append(tuple(map(float, line.strip().split())))
        file.close()
    auc.sort(key=itemgetter(-1))
    order = int(auc[-1][0])
    length = int(auc[-1][1])
    return(length, order)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('train', action='store', help='Path to train data in FASTA format')
    parser.add_argument('output', action='store', help='Output directory to write results')
    parser.add_argument('streme', action='store', help='Directory with STREME motifs. If doesn`t exist it will be created and STREME`s motifs will be calculated')
    parser.add_argument('-t', '--test', action='store', type=str, dest='test',
                        required=False, default='', help='Path to test data in FASTA format. If not specified train data is used.')
    parser.add_argument('-b', '--background', action='store', type=str, dest='background',
                        required=False, default='shuffled', help='Path to background in FASTA format. \
                        It is used for de novo and to estimate pAUC. \
                        if it is not given background is generated by shuffling test data (x5 times)')
    parser.add_argument('-T', '--tmp', action='store', type=str, dest='tmp',
                        required=False, default='./tmp.bamm', help='TMP directory')
    parser.add_argument('-f', '--fpr', action='store', type=float, dest='fpr',
                        required=False, default=0.001, help='The value of FPR to which the pAUC is calculated')
    parser.add_argument('-l', '--min_length', action='store', type=int, dest='min_length',
                        required=False, default=8, help='Minimal length of motif (default value is 8. Don`t use values less than 6)')
    parser.add_argument('-L', '--max_length', action='store', type=int, dest='max_length',
                        required=False, default=20, help='Maximal length of motif (default value is 20. Don`t use values more than 30)')
    parser.add_argument('-s', '--step', action='store', type=int, dest='step',
                        required=False, default=4, help='The step with which the length of the motif will increase (default value is 4)')
    parser.add_argument('-o', '--min_order', action='store', type=int, dest='min_order',
                        required=False, default=1, help='Minimal order of BaMM (default value is 1. Don`t use values less than 1)')
    parser.add_argument('-O', '--max_order', action='store', type=int, dest='max_order',
                        required=False, default=4, help='Maximal order of BaMM (default value is 4. Don`t use values more than 5)')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()

    train_path = args.train
    test_path = args.test
    backgroud_path = args.background
    pwm_dir = args.streme
    tmp_dir = args.tmp
    write_dir = args.output
    fpr_threshold = args.fpr
    min_length = args.min_length
    max_length = args.max_length
    min_order =  args.min_order
    max_order = args.max_order
    step_length = args.step

    if not os.path.exists(train_path):
        print(f'Train file doesn`t exist - {train_path}. Exit')
        sys.exit(1)
    if not step_length >= 0:
        print(f'Warning! Length change step < 1. Actual value is {step_length}. The value will be adjusted to 1.')
        step_length = 1
    if max_length > 30:
        print(f'Warning! Maximal value of length > 40. Actual value is {max_length}. The value will be adjusted to 30.')
        max_length = 30
    if min_length < 6:
        print(f'Warning! Minimal value of length < 6. Actual value is {min_length}. The value will be adjusted to 6.')
        min_length = 6
    if test_path == '':
        test_path = train_path
    learn_optimized_bamm(train_path, test_path, backgroud_path, tmp_dir, write_dir, pwm_dir,
                        fpr_threshold=fpr_threshold, min_length=min_length, max_length=max_length,
                        step_length=step_length, min_order=min_order, max_order=max_order)
    pass


if __name__ == '__main__':
    main()
