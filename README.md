# The boring de novo search for motifs (BDSM)

## Introduction

This repository contains scripts for de novo searching of motifs that were used in the IBIS challange. The repository includes five scripts in total. Two preparatory scripts (1_prepare_data.py, 2_dibackground.py), for de novo search of motifs using the PWM model with STREME[^a] (3_streme.py), for de novo search of motifs using the BaMM[^b] model (4_bamm.py), as well as a technical script (5_scan_best_by_bamm.py), which allows to find the best hits of the BaMM model. The last script is needed to format the results of the AAA model (in our case AAA model is BaMM).

## Requirements

Clone repository:
```console
git clone https://github.com/ubercomrade/bdsm
```

Then create and activate a conda environment:
```console
conda env create -f environment.yml
conda activate bdsm
```

In order to work with the BaMM[^b] model it is necessary to install the corresponding tool: https://github.com/soedinglab/BaMMmotif2. After BaMM is installed, add the path to the executable files to the PATH variable, otherwise script `4_bamm.py` will not work.

## Usage

### Preparation step (1)

_Timing ~1 min_

The command `1_prepare_data.py -h` return:

```
usage: 1_prepare_data.py [-h] [-c COUNT] output peaks genome

positional arguments:
  output                Output directory to write results
  peaks                 Path to PEAKS in narrowPeak format
  genome                Path to GENOME in FASTA format

options:
  -h, --hESDEG includes five commands: `jaspar`, `hocomoco`, `deg`, `set` and `annotaion`.
elp            show this help message and exit
  -c COUNT, --count COUNT
                        The the total number of peaks that will be taken for
                        the test and train samples from the original data.
                        Only even peaks are taken according to their order for
                        test sample. Only odd peaks are taken according to
                        their order for train sample.
```

This script creates a directory and subdirectories necessary to run other programs, and also divides the initial peak sampling into tour sampling and training sampling.

**First positional argument** `output`:

Name of the directory where all results and other data will be put

**Second positional argument** `peaks`:

Path to peaks in narrowPeak format (see MACS2 or UCSC)

**Third positional argument** `genome`:

Path to reference genome in FASTA format

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-c; --COUNT` :

The the total number of peaks that will be taken for the test and train samples from the original data. Only even peaks are taken according to their order for test sample. Only odd peaks are taken according to their order for train sample.

### Preparation (2)

Generating dinucleotide background. To speed up this step PyPy can be used instead of CPython.

_Timing 10–120 min and more_

The command `2_dibackground.py -h` return:

```
usage: 2_dibackground.py [-h] [-thr THR_CONTENT] [-c COUNTS]
                         output fasta genome

positional arguments:
  output                Path to write output file in FASTA format
  fasta                 Path to FASTA with head structure: >CHR:START-END
  genome                Path to GENOME in FASTA format

options:
  -h, --help            show this help message and exit
  -thr THR_CONTENT, -thr_content THR_CONTENT
                        The value of maximum difference in GC composition
                        between the original peak and the generated peak.
                        Default value is 0.01
  -c COUNTS, --counts COUNTS
                        Number of peaks that will be generated for one peak
                        from FASTA. Default value is 5
```

#### Required arguments description

**First positional argument** `output`:

The path of background to be recorded in FASTA format

**Second positional argument** `fasta`:

The path to the tutorial/test sample in FASTA format. They are created with the first command (1_prepare_data.py), so see the results of the first command.

**Third positional argument** `genome`:

Path to reference genome in FASTA format

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-thr` :

The value of maximum difference in GC composition between the original peak and the generated peak. Default value is 0.01

**Third optional argument** `-c; --COUNT` :

Number of peaks that will be generated for one peak from FASTA (train/test sample). Default value is 5

### de novo PWM (3)

Running STREME[^a] to search for de novo motifs.

The command `3_streme.py -h` return:

```
usage: 3_streme.py [-h] [-t TEST] [-b BACKGROUND] [-T TMP] [-f FPR]
                   [-l MIN_LENGTH] [-L MAX_LENGTH] [-s STEP]
                   train output

positional arguments:
  train                 Path to train data in FASTA format
  output                Output directory to write results

options:
  -h, --help            show this help message and exit
  -t TEST, --test TEST  Path to test data in FASTA format. If not specified
                        train data is used.
  -b BACKGROUND, --background BACKGROUND
                        Path to background in FASTA format. It is used for de
                        novo and to estimate pAUC. if it is not given
                        background is generated by shuffling test data (x5
                        times)
  -T TMP, --tmp TMP     TMP directory
  -l MIN_LENGTH, --min_length MIN_LENGTH
                        Minimal length of motif (default value is 8. Don`t use
                        values less than 6)
  -L MAX_LENGTH, --max_length MAX_LENGTH
                        Maximal length of motif (default value is 20. Don`t
                        use values more than 30)
  -s STEP, --step STEP  The step with which the length of the motif will
                        increase (default value is 4)
```

#### Required arguments description

**First positional argument** `train`:

The path of train sample in FASTA format (see 1_prepare_data.py output)

**Second positional argument** `output`:

The path to write results (directory).

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-t, --test` :

The path of test sample in FASTA format (see 1_prepare_data.py output)

**Third optional argument** `-b; --background` :

The path of background in FASTA format (see 2_dibackground.py output)

**Fourth optional argument** `-T; --tmp` :

The path of TMP directory

**Fifth optional argument** `-l; --min_length` :

Minimal length of motif

**Sixth optional argument** `-L; --max_length` :

Maximal length of motif

**Seventh optional argument** `-L; --max_length` :

The step with which the length of the motif will increase


### de novo BaMM (4)

Running BaMM[^b] to search for de novo motifs. To speed up this step PyPy can be used instead of CPython.

The command `4_bamm.py -h` return:

```
usage: 4_bamm.py [-h] [-t TEST] [-b BACKGROUND] [-T TMP] [-f FPR]
                 [-l MIN_LENGTH] [-L MAX_LENGTH] [-s STEP] [-o MIN_ORDER]
                 [-O MAX_ORDER]
                 train output streme

positional arguments:
  train                 Path to train data in FASTA format
  output                Output directory to write results
  streme                Directory with STREME motifs. If doesn`t exist it will
                        be created and STREME`s motifs will be calculated

options:
  -h, --help            show this help message and exit
  -t TEST, --test TEST  Path to test data in FASTA format. If not specified
                        train data is used.
  -b BACKGROUND, --background BACKGROUND
                        Path to background in FASTA format. It is used for de
                        novo and to estimate pAUC. if it is not given
                        background is generated by shuffling test data (x5
                        times)
  -T TMP, --tmp TMP     TMP directory
  -f FPR, --fpr FPR     The value of FPR to which the pAUC is calculated
  -l MIN_LENGTH, --min_length MIN_LENGTH
                        Minimal length of motif (default value is 8. Don`t use
                        values less than 6)
  -L MAX_LENGTH, --max_length MAX_LENGTH
                        Maximal length of motif (default value is 20. Don`t
                        use values more than 30)
  -s STEP, --step STEP  The step with which the length of the motif will
                        increase (default value is 4)
  -o MIN_ORDER, --min_order MIN_ORDER
                        Minimal order of BaMM (default value is 1. Don`t use
                        values less than 1)
  -O MAX_ORDER, --max_order MAX_ORDER
                        Maximal order of BaMM (default value is 4. Don`t use
                        values more than 5)
```

#### Required arguments description

**First positional argument** `train`:

The path of train sample in FASTA format (see 1_prepare_data.py output)

**Second positional argument** `output`:

The path to write results (directory).

**Third positional argument** `streme`:

The path with streme results (directory created by 3_streme.py).

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT

**Second optional argument** `-t, --test` :

The path of test sample in FASTA format (see 1_prepare_data.py output)

**Third optional argument** `-b; --background` :

The path of background in FASTA format (see 2_dibackground.py output)

**Fourth optional argument** `-T; --tmp` :

The path of TMP directory

**Fifth optional argument** `-l; --min_length` :

Minimal length of motif

**Sixth optional argument** `-L; --max_length` :

Maximal length of motif

**Seventh optional argument** `-L; --max_length` :

The step with which the length of the motif will increase

**Eighth optional argument** `-o; --min_order` :

Minimal value of order of BaMM model

**Ninth optional argument** `-O; --max_order` :

Maximal value of order of BaMM model

### de novo BaMM (4)

Scan FASTA by BaMM[^b] to find best hits in each sequences. To speed up this step PyPy can be used instead of CPython.

The command `5_scan_best_by_bamm.py -h` return:

```
usage: 5_scan_best_by_bamm.py [-h] bamm fasta output tf

positional arguments:
  bamm        Dir with BaMM model
  fasta       path to FASTA
  output      Path to write results
  tf          Name of tf

options:
  -h, --help  show this help message and exit
```

#### Required arguments description

**First positional argument** `bamm`:

The directory with BaMM model (see 4_bamm.py output)

**Second positional argument** `fasta`:

The path of FASTA file to scan.

**Third positional argument** `output`:

The path to write results in IBIS format.

Example of output:
```
tag	LEUTX
lzxxukv	0.70490
oynfkok	0.87811
lyswvww	0.73776
rxfymzl	0.73776
zojmjmt	0.70774
emvolbt	0.76957
vbfoatm	0.76957
djaudvh	0.76957
acrewpu	0.76957
```

**Fourth positional argument** `streme`:

The name of TF. It will be used as name of column (see output example)

#### Optional arguments description

**First optional argument** `-h; --help` :

Print help to STDOUT


[^a]: Bailey, T. L. (2021). STREME: accurate and versatile sequence motif discovery. Bioinformatics, 37(18), 2834–2840. https://doi.org/10.1093/bioinformatics/btab203
[^b]: Ge, W., Meier, M., Roth, C., & Söding, J. (2021). Bayesian Markov models improve the prediction of binding motifs beyond first order. NAR Genomics and Bioinformatics, 3(2), lqab026. https://doi.org/10.1093/nargab/lqab026
