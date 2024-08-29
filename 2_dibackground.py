#!/usr/bin/env python3

import sys
import argparse
import random
import itertools
from pyfaidx import Fasta, Faidx
from math import sqrt, dist

def read_fasta(path):
    fasta = []
    with open(path) as file:
        for line in file:
            if line.startswith('>'):
                chromosome, start_end = line[1:].strip().split(':')
                start, end = map(int, start_end.split('-'))
            else:
                fasta.append(((chromosome, start, end), line.strip().upper()))
    return fasta


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def wrtie_background(path, background):
    with open(path, 'w') as file:
        for info, seq in background:
            chromosome, start, end = info
            file.write(f'>{chromosome}:{start}-{end}\n')
            file.write(f'{seq}\n')
    pass


# def calculate_dicontent(seq, dinucleotides):
#     length = len(seq)
#     content = {di:0 for di in dinucleotides}
#     for i in range(length - 1):
#         di = seq[i:i+2]
#         if di in content:
#             content[di] += 1
#     for di in content.keys():
#         content[di] = content[di] / length
#     return content


def calculate_dicontent(seq, dinucleotides):
    length = len(seq) - seq.count('N')
    return [seq.count(di) / length for di in dinucleotides]


def calculate_rmsd(vec1, vec2):
    length = len(vec1)
    return sqrt(sum([(i - j)**2 for i,j in zip(vec1, vec2)]) / length)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('output', action='store', help='Path to wrtite output file in FASTA format')
    parser.add_argument('fasta', action='store', help='Path to FASTA with head strucutre: >CHR:START-END')
    parser.add_argument('genome', action='store', help='Path to GENOME in FASTA format')
    parser.add_argument('-thr', '-thr_content', action='store', type=float, dest='thr_content',
                        required=False, default=0.01, help='The value of maximum difference in GC composition \
                        between the original peak and the generated peak. Default value is 0.01')
    parser.add_argument('-c', '--counts', action='store', type=int, dest='counts',
                        required=False, default=5, help='Number of peaks that will be generated for one peak from FASTA. Default value is 5')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()

    output = args.output
    genome_path = args.genome
    fasta_path = args.fasta
    thr_content = args.thr_content
    counts = args.counts

    genome = Fasta(genome_path)
    idx = Faidx(genome_path)
    chromosomes = list(idx.index.keys())
    chromosomes = [i for i in chromosomes if 'chr' in i]

    genome_mem = dict()
    for chr in chromosomes:
        genome_mem[chr] = genome[chr][:]

    fasta = read_fasta(fasta_path)

    background = []
    dinucleotides = [''.join(i) for i in itertools.product('ACGT', repeat=2)]
    number_of_dinucleotides = len(dinucleotides)
    for info, seq in fasta:
        length = len(seq)
        if seq.count('N') == length:
            continue

        dicontent = calculate_dicontent(seq, dinucleotides)
        #print(dicontent)
        sbackground = []
        for i in range(500000):
            chromosome = random.choice(chromosomes)
            if idx.index[chromosome]['rlen'] > length * 2:
                start = random.randint(0, idx.index[chromosome]['rlen'] - length)
                end = start + length
                if chromosome == info[0]:
                    if info[2] > start and info[1] < end:
                        continue
                rnd_seq = genome_mem[chromosome][start:end]
                rnd_seq = rnd_seq.seq.upper()
                if rnd_seq.count('N') / len(rnd_seq) > 0.2: continue


                #approach 1
                # rnd_dicontent = calculate_dicontent(rnd_seq, dinucleotides)
                # passed = len(list(itertools.takewhile(lambda di: abs(1 - rnd_dicontent[di] / dicontent [di]) < thr_content, dinucleotides)))
                # rnd_dicontent = calculate_dicontent(complement(rnd_seq), dinucleotides)
                # passed_rc = len(list(itertools.takewhile(lambda di: abs(1 - rnd_dicontent[di] / dicontent [di]) < thr_content, dinucleotides)))

                #approach 2
                rnd_dicontent = calculate_dicontent(rnd_seq, dinucleotides)
                #distance = calculate_rmsd(dicontent, rnd_dicontent)
                distance = dist(dicontent, rnd_dicontent)
                rnd_dicontent = calculate_dicontent(complement(rnd_seq), dinucleotides)
                #distance_rc = calculate_rmsd(dicontent, rnd_dicontent)
                distance_rc = dist(dicontent, rnd_dicontent)

                #if passed == number_of_dinucleotides or passed_rc == number_of_dinucleotides:
                if distance < thr_content or distance_rc < thr_content:
                    sbackground.append(((chromosome, start, end), rnd_seq))
                    if len(sbackground) == counts:
                        break
        background += sbackground
    print(len(background))
    wrtie_background(output, background)
    pass


if __name__ == '__main__':
    main()
