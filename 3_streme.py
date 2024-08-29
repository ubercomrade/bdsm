#!/usr/bin/env python3

import os
import sys
import random
import shutil
import bisect
import argparse
import subprocess
import operator
import numpy as np
from shutil import copyfile
from operator import itemgetter
from numpy.lib.stride_tricks import sliding_window_view


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


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


def seq_to_int(sequenes):
    converter = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    results = []
    for seq in sequenes:
        vector = np.array([converter[i] for i in seq])
        results.append(vector)
    return(results)


def dict_to_array(motif):
    motif = [motif[i] for i in motif.keys()]
    return np.array(motif)


def array_to_dict(array):
    nucleotides = ['A', 'C', 'G', 'T']
    motif = {nuc:list(array[index]) for index, nuc in enumerate(nucleotides)}
    return motif


def pfm_to_pwm(pfm):
    background = 0.25
    pwm = np.log2(pfm / background)
    return pwm


# def pfm_to_pwm_cisbp(pfm):
#     background = 0.25
#     pwm = np.log2((pfm + 0.01) / background)
#     return pwm


# def pcm_to_pfm(pcm):
#     number_of_sites = pcm.sum(axis=0)
#     nuc_pseudo = 0.25
#     pfm = (pcm + nuc_pseudo) / (number_of_sites + 1)
#     return pfm


def min_score(matrix):
    return np.sum(matrix.min(axis=0))


def max_score(matrix):
    return np.sum(matrix.max(axis=0))


def to_score(matrix, norm_value):
    min_s = min_score(matrix)
    max_s = max_score(matrix)
    score = norm_value * (max_s - min_s) + min_s
    return score


def scores_of_sites(pwm, matrix_of_sites):
    return np.sum(np.take_along_axis(pwm, matrix_of_sites, axis=0), axis=1)


def streme_to_pwm(pfm):
    pfm = dict_to_array(pfm)
    pwm = pfm_to_pwm(pfm)
    return pwm


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


def all_scores_pwm(numeric_sequenes, pwm, length_of_site):
    # scores, scores_rc = []
    # append, append_rc = scores.append, scores_rc.append
    scores = []
    append = scores.append
    for i in range(len(numeric_sequenes)):
        matrix_of_sites = sliding_window_view(numeric_sequenes[i], length_of_site)
        matrix_of_sites_rc = np.flip(np.array((matrix_of_sites -3) * -1), axis=1)
        append(scores_of_sites(pwm, matrix_of_sites))
        append(scores_of_sites(pwm, matrix_of_sites_rc))
    return np.concatenate(scores)


def best_scores_pwm(numeric_sequenes, pwm, length_of_site):
    # scores, scores_rc = []
    # append, append_rc = scores.append, scores_rc.append
    scores = []
    append = scores.append
    for i in range(len(numeric_sequenes)):
        matrix_of_sites = sliding_window_view(numeric_sequenes[i], length_of_site)
        matrix_of_sites_rc = np.array((matrix_of_sites -3) * -1)[::-1]
        score = np.max(scores_of_sites(pwm, matrix_of_sites))
        score_rc = np.max(scores_of_sites(pwm, matrix_of_sites_rc))
        scores.append(np.max([score, score_rc]))
    return scores


def write_roc(path, data):
    header = list(data.keys())
    with open(path, 'w') as file:
        file.write('\t'.join(header) + '\n')
        for i in zip(*data.values()):
            file.write('\t'.join((map(str,i))) + '\n')
    return(0)


def write_table(path, *args):
    with open(path, 'a') as file:
        file.write('\t'.join(list(map(str, args))[::-1])+'\n')
    return(0)


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


# def calculate_roc(fprs):
#     table = {'TPR': [0], 'FPR': [0], 'SITES': [0]}
#     current_number_of_sites = 0
#     total_number_of_sites = len(fprs)
#     for i in range(0, total_number_of_sites):
#         if fprs[i] > fprs[i - 1]:
#             table['TPR'].append(current_number_of_sites / total_number_of_sites)
#             table['FPR'].append(fprs[i - 1])
#             table['SITES'].append(current_number_of_sites)
#         current_number_of_sites += 1
#     table['TPR'].append(current_number_of_sites / total_number_of_sites)
#     table['FPR'].append(fprs[i])
#     table['SITES'].append(current_number_of_sites)
#     return(table)
#
#
# def calculate_roc(true_scores, false_scores):
#     tprs = []
#     fprs = []
#     true_scores.sort()
#     false_scores.sort()
#     true_scores_uniq = list(set(true_scores))
#     true_scores_uniq.sort(reverse=True)
#     false_length = len(false_scores)
#     true_length = len(true_scores)
#     for score in true_scores_uniq:
#         tpr = (true_length - bisect.bisect_right(true_scores, score)) / true_length
#         fpr = (false_length - bisect.bisect_right(false_scores, score)) / false_length
#         if fpr == 0:
#             fpr = 0.5 / false_length
#         tprs.append(tpr)
#         fprs.append(fpr)
#     return(tprs, fprs)
#
#
# def creat_table_bootstrap(true_scores, false_scores):
#     cdef list table = []
#     true_scores.sort(reverse=True)
#     false_scores.sort(reverse=True)
#     false_length = len(false_scores)
#     true_length = len(true_scores)
#     for tpr in [round(i * 0.01, 2) for i in range(5,105, 5)]:
#         score = true_scores[round(true_length * tpr) - 1]
#         actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
#         fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
#         table.append({'Scores': score, 'TPR': tpr, 'ACTUAL_TPR': actual_tpr, 'FPR': fpr})
#     return(table)


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


def run_streme_hmm_background(fasta_path, dir_out, motif_length):
    args = ['streme', '--p', fasta_path,
            '--kmer', '4',
           '--oc', dir_out,
           '--objfun', 'de',
           '--w', str(motif_length),
           '-nmotifs', '5']
    p = subprocess.run(args, shell=False, capture_output=True)
    return(0)


def run_streme(fasta_path, backgroud_path, dir_out, motif_length):
    args = ['streme', '--p', fasta_path,
            '--n', backgroud_path,
           '--oc', dir_out,
           '--objfun', 'de',
           '--w', str(motif_length),
           '--nmotifs', '5',
           '--neval', '20',
           '--nref', '4',
           '--niter', '20']
    p = subprocess.run(args, shell=False, capture_output=True)
    return(0)


# def run_meme(fasta_path, backgroud_path, dir_out, motif_length):
#     args = ['streme', fasta_path,
#             '-neg', backgroud_path,
#            '--oc', dir_out,
#            '-objfun', 'de',
#            '-w', str(motif_length),
#            '-mod', 'zoops',
#            '-nmotifs', '5']
#     p = subprocess.run(args, shell=False, capture_output=False)
#     return(0)


def run_streme_hmm_background(fasta_path, dir_out, motif_length):
    args = ['streme', '--p', fasta_path,
            '--order', '5',
           '--oc', dir_out,
           '--objfun', 'de',
           '--w', str(motif_length),
           '-nmotifs', '5']
    p = subprocess.run(args, shell=False, capture_output=True)
    return(0)


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


def write_pwm(output, tag, pwm):
    with open(output + '/' + tag + '.pwm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(pwm['A'], pwm['C'], pwm['G'], pwm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))


def write_pfm(output, tag, pfm):
    with open(output + '/' + tag + '.pfm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\n'.format(i[0], i[1], i[2], i[3]))


def learn_optimized_pwm(train_path, test_path, backgroud_path, tmp_dir, write_dir,
                        fpr_threshold=0.001, min_length=8, max_length=20, step_length=4):

    # Read peaks and background
    # train_numeric = seq_to_int(read_fasta(train_path))
    # test_numeric = seq_to_int(read_fasta(test_path))
    if os.path.isfile(backgroud_path):
        test_numeric = seq_to_int(read_fasta(test_path))
        background_numeric = seq_to_int(read_fasta(backgroud_path))
    else:
        test_numeric = seq_to_int(read_fasta(test_path))
        indexes = np.random.choice(np.arange(0, len(test_numeric)), len(test_numeric) * 5)
        background_numeric = [np.copy(test_numeric[i]) for i in indexes]
        for i in background_numeric:
            np.random.shuffle(i)

    for length in range(min_length, max_length + 1, step_length):
        #Run STREME
        if os.path.isfile(backgroud_path):
            run_streme(train_path, backgroud_path, tmp_dir, length)
        else:
            run_streme_hmm_background(train_path, tmp_dir, length)

        #Read STREME motif
        pfm, background, length, nsites = parse_streme(tmp_dir + '/streme.txt')
        pwm = streme_to_pwm(pfm)
        pwm = np.concatenate((pwm, np.array([pwm.shape[1] * [-100]])), axis=0)

        #Write motifs
        tag = f'motif_{length}'
        write_meme(f'{write_dir}/statistics/', tag, pfm, background, nsites)
        write_pwm(f'{write_dir}/statistics/', tag, array_to_dict(pwm))
        write_pfm(f'{write_dir}/statistics/', tag, pfm)
        copyfile(f'{tmp_dir}/streme.html', f'{write_dir}/statistics/{tag}.html')

        # vs frequencies of peaks
        true_scores = best_scores_pwm(test_numeric, pwm, length)
        false_scores = best_scores_pwm(background_numeric, pwm, length)

        #Calculate FPRs -> ROC -> AUC
        fprs = calculate_fprs(true_scores, false_scores)
        merged_roc = calculate_merged_roc(fprs)
        auc = calculate_particial_auc(merged_roc['TPR'], merged_roc['FPR'], 1.01)

        prc = calculate_prc(true_scores, false_scores)
        auc_prc = calculate_particial_auc(prc['PRECISION'], prc['RECALL'], 1.01)

        # vs frequencies of sites
        #Scan peaks by motif
        true_scores = best_scores_pwm(test_numeric, pwm, length)
        false_scores = all_scores_pwm(background_numeric, pwm, length)

        #Calculate FPRs -> ROC -> pAUC
        fprs = calculate_fprs(true_scores, false_scores)
        staged_roc = calculate_staged_roc(fprs, step=1)
        merged_roc = calculate_merged_roc(fprs)
        pauc = calculate_particial_auc(merged_roc['TPR'], merged_roc['FPR'], fpr_threshold)

        #Write ROC and pAUC
        print(f"Length {length}; pAUC at {fpr_threshold} = {pauc}; aucROC = {auc:.4f}; aucPRC = {auc_prc:.4f}")
        write_table(f'{write_dir}/statistics/auc.txt', (pauc / 0.001 + auc + auc_prc), length)
        write_roc(f'{write_dir}/statistics/roc_merged_{length}.txt', merged_roc)
        write_roc(f'{write_dir}/statistics/roc_{length}.txt', staged_roc)
    shutil.rmtree(tmp_dir)
    return 0


def choose_best_model(path):
    auc = []
    with open(path) as file:
        for line in file:
            auc.append(tuple(map(float, line.strip().split())))
        file.close()
    auc.sort(key=itemgetter(-1))
    length = int(auc[-1][0])
    return(length)


def de_novo_with_oprimization_pwm_streme(train_path, test_path, backgroud_path, tmp_dir, write_dir,
                                         fpr_threshold=0.001, min_length=8, max_length=20, step_length=4):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    if not os.path.isdir(write_dir):
        os.mkdir(write_dir)
    if not os.path.isdir(f'{write_dir}/statistics'):
        os.mkdir(f'{write_dir}/statistics')
    if os.path.exists(f'{write_dir}/statistics/auc.txt'):
        os.remove(f'{write_dir}/statistics/auc.txt')
    learn_optimized_pwm(train_path, test_path, backgroud_path, tmp_dir, write_dir,
                        fpr_threshold, min_length, max_length, step_length)
    length = choose_best_model(f'{write_dir}/statistics/auc.txt')
    print(f'Model have best perfomance with length = {length}')
    copyfile(f'{write_dir}/statistics/roc_{length}.txt',
             f'{write_dir}/roc.txt')
    copyfile(f'{write_dir}/statistics/roc_merged_{length}.txt',
             f'{write_dir}/roc_merged_{length}.txt')
    copyfile(f'{write_dir}/statistics/motif_{length}.pfm',
             f'{write_dir}/motif.pfm')
    copyfile(f'{write_dir}/statistics/motif_{length}.pwm',
             f'{write_dir}/motif.pwm')
    copyfile(f'{write_dir}/statistics/motif_{length}.meme',
             f'{write_dir}/motif.meme')
    return 0


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('train', action='store', help='Path to train data in FASTA format')
    parser.add_argument('output', action='store', help='Output directory to write results')
    parser.add_argument('-t', '--test', action='store', type=str, dest='test',
                        required=False, default='', help='Path to test data in FASTA format. If not specified train data is used.')
    parser.add_argument('-b', '--background', action='store', type=str, dest='background',
                        required=False, default='shuffled', help='Path to background in FASTA format. \
                        It is used for de novo and to estimate pAUC. \
                        if it is not given background is generated by shuffling test data (x5 times)')
    parser.add_argument('-T', '--tmp', action='store', type=str, dest='tmp',
                        required=False, default='./tmp.streme', help='TMP directory')
    # parser.add_argument('-f', '--fpr', action='store', type=float, dest='fpr',
    #                     required=False, default=0.001, help='The value of FPR to which the pAUC is calculated')
    parser.add_argument('-l', '--min_length', action='store', type=int, dest='min_length',
                        required=False, default=8, help='Minimal length of motif (default value is 8. Don`t use values less than 6)')
    parser.add_argument('-L', '--max_length', action='store', type=int, dest='max_length',
                        required=False, default=20, help='Maximal length of motif (default value is 20. Don`t use values more than 30)')
    parser.add_argument('-s', '--step', action='store', type=int, dest='step',
                        required=False, default=4, help='The step with which the length of the motif will increase (default value is 4)')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()

    train_path = args.train
    test_path = args.test
    backgroud_path = args.background
    tmp_dir = args.tmp
    write_dir = args.output
    #fpr_threshold = args.fpr
    min_length = args.min_length
    max_length = args.max_length
    step_length = args.step

    fpr_threshold = 0.001
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

    de_novo_with_oprimization_pwm_streme(train_path, test_path, backgroud_path, tmp_dir, write_dir,
                                     fpr_threshold=fpr_threshold, min_length=min_length, max_length=max_length, step_length=step_length)
    pass


if __name__ == '__main__':
    main()
