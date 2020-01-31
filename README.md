[![travis](https://img.shields.io/travis/kimlaborg/kmtools.svg?style=flat-square)](https://travis-ci.org/kimlaborg/NGSKit/)
[![codecov](https://img.shields.io/codecov/c/github/kimlaborg/kmtools.svg?style=flat-square)](https://codecov.io/gh/kimlaborg/NGSKit)

# NGSKit

This repository contains the full collection of tools of our lab preprocessing NGS data pipeline. These tools are heavily custom to Kimlab's screenings, consequently please take a look at the references section (above) to read in more detail about some application cases of this package.

The package can be divided into three different sets of tools; NGS preprocessing data, secondly tools and pipelines for analysis of the screening, and finally, the library encoding tools.

## NGSkit utilized on the subsequent papers.

**A PxL motif promotes timely cell cycle substrate dephosphorylation by the Cdc14 phosphatase.** Kataria M, Mouilleron S, Seo MH, Corbi-Verge C, Kim PM, Uhlmann F. Nat Struct Mol Biol. 2018 Nov 19. doi: 10.1038/s41594-018-0152-3. [DATA](https://zenodo.org/record/3633357#.XjSXkuF7lGo)

**A multi-reporter bacterial 2-hybrid assay for the fast, simple, and dynamic assay of PDZ domain – peptide interactions.**
Ichikawa, David; Corbi-Verge, Carles; Shen, Michael; Snider, Jamie; Wong, Victoria; Stagljar, Igor; Kim, Philip; Noyes, Marcus. ACS Synth Biol. 2019 Apr 18. doi: 10.1021/acssynbio.8b00499.
[DATA](https://zenodo.org/record/2580337#.Xff0VtF7mV4)


# Installing

```bash
pip install -r requirements.txt   . 
```

## Requirements

  - pandas
  - numpy
  - biopython
  - scipy
  - pytest
  - scikit-learn
  - tqdm
  - python-Levenshtein
  - weblogo
  - matplotlib


## Optionals

  - Pandaseq
  - Fastqc



# NGS Preprocessing Data.

Screenings samples are usually pooled sequenced using illimuna sequencing technology. Each sample is barcoded with a short nucleotide sequence. The preprocessing is divided into two independent steps, demultiplexing and trimming. This design was chossed to facilitate troubleshooting.  The first step is to read the output file from the sequencing filtering low-quality reads. Files can be demultiplexed on the spot, but it can help full under certain circumstances to have more control of this step. 
During the Trimming step, the target sequence is extracted from the flanking sequences to troubleshoot. A few extra options can be activate,  allowing, for instance, dynamic length extraction. 


##  Demultiplexation

Demultiplexation Fastq sequences tool:

```bash
demultiplexer  -b [BarCode_file.inp] -i [deep_seq_file.fastq] -o [folder_name] -l 54 -m QUICK --misreads_cutoff_cons 2


optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FASTQS [INPUT_FASTQS ...], --input_fastqs INPUT_FASTQS [INPUT_FASTQS ...]
                        input_fastqs FASTQ file or files (demultiplex)
  -b BARCODE_FILE, --barcode_file BARCODE_FILE
                        File that contains barcodes and cosntant regions
  -o OUT_DIR, --out_dir OUT_DIR
                        Output folder, called demultiplex by default
  -m {quick,standard,simple,dynamic}, --demultiplexation_method {quick,standard,simple,dynamic}
                        Type of demultiplexation by default; STANDARD `quick`:
                        Only the first barcode and constant region will be
                        check `standard`: Both barcodes and constant regions
                        will be check `simple`: Only the barcodes are used
                        `dynamic`: frame shift search, Flexible search of the
                        second constant region and barcode
  --misreads_cutoff_cons MISREADS_CUTOFF_CONS
                        Max number of misreading allowed in the constant
                        constant_region (default 2)
  --misreads_cutoff_barcode MISREADS_CUTOFF_BARCODE
                        Max number of misreading allowed in the constant
                        constant_region (default 1)
  --dump                Dump constant regions
  --no-save_frequencies
                        Do not Save match frequencies
```

## Trimming

Trimming Fastq sequences tool Usage Trimming:

```bash
trimmer  -d [demultiplexedFolder]-b [BarCode_file.inp] -q [Quality threshold] -m [method] --output_fmt fasta

optional arguments:
  -h, --help            show this help message and exit
  -d INPUT_FOLDER, --input_folder INPUT_FOLDER
                        Folder contains demultiplexed folders and files
  -b BARCODE_FILE, --barcode_file BARCODE_FILE
                        File that contains barcodes and cosntant regions
  -o OUT_FOLDER, --out_folder OUT_FOLDER
                        Output folder, called Sequences by default
  -m {standard,dynamic}, --trimming_method {standard,dynamic}
                        standard Trimm sequences according barcode file
                        configuration, ignores float window output files
                        dynamic Trimm sequences using file lenght label, or
                        output of float window demultiplex
  -q QUALITY, --quality QUALITY
                        Quality reading threshold (default 30)
  --output_fmt OUTPUT_FMT
                        Output format, default fasta
  --force-lenght FORCE_LENGHT
                        force a lenght and ignore file label, overwrites
                        dynamic option

```

## Barcodes

Small utility to setup sample barcodes into the format to be used by the demultiplexer tool. Also, it has a few small functionalities to modify the sequences to transform sequences to reverse and complementary, etc. 

```bash

usage: barcodes [-h] -b BARCODE_FILE [-o OUT_PREFIX] [--rc_b2] [--rc_c2]
                [--skip_header]

<<<<<<< HEAD
Barcodes tool Usage : %prog -b [BarCode_file.excel] -o [to_demultiplex_]

optional arguments:
  -h, --help            show this help message and exit
  -b BARCODE_FILE, --barcode_file BARCODE_FILE
                        File that contains barcodes and cosntant regions
  -o OUT_PREFIX, --out_prefix OUT_PREFIX
                        Output prefix name to_demultiplex by default
  --rc_b2               Reverse and Complemete barcode 2
  --rc_c2               Reverse and Complemete Constant region 2
  --skip_header         Input file has a header, skip it
```


# Library generator

This module contains functions and tools to translate a library of peptides or protein into a nucleotide library. Among other options, the library can be encoded with custom flanking adaptors, restriction enzymes cleaving sites and using preferred codon usages.Check the notebooks folder for some examples. 

# Analysis Tools

=======
>>>>>>> e0af60ecdde084327d3d9a1a04c630cfc4419b6f

# Acknowledgement

This code have been inspired by the combination of script wrote by Juhyun Jeon. I incorporate also a lot of valuable feedback from Satra Nim and Alexey Strokach. Everything under the supervision of Prof. Philip M Kim. 
