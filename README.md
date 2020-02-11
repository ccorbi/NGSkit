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


# Install

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



# NGSkit preprocessing data tools.

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

  -m, --demultiplexation_method {quick,standard,simple,dynamic}
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

  --late_end LATE_END   Number of positions to keep looking for the 2nd
                        constant region after the target lenght. Useful to
                        skip stop codons, or insertions in the designed Seq.
                        Only applys on dynamic method

  --early_end EARLY_END
                        Number of positions to start looking for the 2nd
                        constant. Uset when there are deletions in the
                        designed seq. Only applys on dynamic method

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


# Analysis 

- **TBC**

# Preprocessing Example

**Barcodes and QC**

Before running the entire pipeline, it is strongly recommended to check the quality of the raw data. It is planned to incorporate this functionality in demultiplexing utility, but in the meantime, you can use a third-party tool. There are many tools to get a report on your file. For instance [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
Also, determine is you need to stitch both reverse and forward Fastq files. This is not yet implemented, but again you can use a tool like [Pandaseq](https://github.com/neufeld/pandaseq). 

- *Example pandaseq*

```bash

pandaseq  -f lab_S1_R1_001.fastq.gz -r lab_S1_R2_001.fastq.gz -A ea_util  -F -w merged_file_20180303.fastq -l merge.log 
```

 - *Example barcode file*
  
Before demultiplexing, we need to create a text file with the information of each barcode. 


|Sample name   | Barcode 5'   |Constant Region 5'   | Constant Region 3'  | Barcode 3'   | Design region length   |
|---|---|---|---|---|---|
| Sample_1  | CATA  | TGTGGGTC   | AGGATG  | CATT  | 21   |
| Sample_2  |  CTTT  | TGTGGGTC   | AGGATG  | CATT  | 21 |


Sample name will be use to identify files and folders,so need to be unique  and spaces in the name are not allowed. In the case, there is not a secondary barcode or constant regions; it is needed to add a dash '-' to the empty field. All the barcodes need to be in 5` -> 3'.

```bash
barcodes -b PDZ_illumina.xlsx -o test_PDZ --skip_header
```

The easiest way to do it is to create an excel file with all barcodes with the latter format and run the barcode utility to generate the text file. By default, the utility will make an individual file for a sample. This behaviour is intended to embarrassing parallelize the demultiplexing process by running an independent process over the same raw data but can be modified. If you need another approach, please check the input parameters of both demultiplex and barcode utilities. 


**Running multiple demultiplex process and Trimming**

The idea is to submit multiple jobs, each with an individual sample, in an HPC environment with SLURM.


```bash
for i in *.barcode;do sbatch --partiotion=cpu --mem=4G --wrap="demultiplexer -i Lane1_S1_L001_R1_001.fastq.gz  Lane2_S2_L002_R1_001.fastq.gz Lane3_S3_L003_R1_001.fastq.gz Lane4_S4_L004_R1_001.fastq.gz -b ${i} -m quick -o my_folder_demultiplex_data";done
```

If you run the pipeline in a local machine with multiple cores (i.e 8), you can run multiple instances of the demultiplexer too. 

```bash
parallel -j 8  demultiplexer -i Lane1_S1_L001_R1_001.fastq.gz  Lane2_S2_L002_R1_001.fastq.gz Lane3_S3_L003_R1_001.fastq.gz Lane4_S4_L004_R1_001.fastq.gz -b {} -m quick -o my_folder_demultiplex_data ::: *.barcode
```

The output of this step is a folder with the name of the sample and within one or multiple Fastq files.  The number on the file shows the length of the design region, this is important when there are multiple design regions, or we allow for insertion and deletions in the design region and run demultiplexing with the dynamic method. 

There is a folder named Stats that contains counting of sequences for each level of match, barcode, barcode and the constant region, design region length, etc. This log can help to detect insertions and deletions in barcodes, or other errors. I recommend checking the file after the job is done. Next, we need to trim the design regions from the demultiplexed sequences. This extra step is a separate of demultiplexing to help to debug potential problems and avoid to demultiplex multiple times with different paraments.  However, everything can be packed in a single job. 

```bash
for i in *.barcode;do trimer -d ./my_folder_demultiplex_data/ -b ${i}; done
```

Finally, depending on the type of screening and analysis, you want to do different counting, take a look at pipelines and counter_reads_entropy.py  and counter_reads.py . 

**Mapping to a library**

There are multiple approaches to remove sequences not present in the designed library. One approach is to use [Bowtie](http://bowtie-bio.sourceforge.net/manual.shtml).  Set up the parameters according to your needs. 


```bash
bowtie --best -v 0 ../../Libs/disdlib_4_endomembrane_system -f Disorderome_lib_naive_4.fasta > ./bowtie/Disorderome_lib_naive_4.out
```

You can parse bowtie output using counter_map_reads.py  or any other custom tool. 



# Acknowledgements

This code have been inspired by the combination of script wrote by Juhyun Jeon. I incorporate also a lot of valuable feedback from Satra Nim and Alexey Strokach. Everything under the supervision of Prof. Philip M Kim. 
