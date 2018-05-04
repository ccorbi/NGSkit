
[![travis](https://img.shields.io/travis/kimlaborg/kmtools.svg?style=flat-square)](https://travis-ci.org/kimlaborg/NGSKit/)
[![codecov](https://img.shields.io/codecov/c/github/kimlaborg/kmtools.svg?style=flat-square)](https://codecov.io/gh/kimlaborg/NGSKit)

# NGSKit

Collection basic tools to preprocess NGS data for the lab Pipelines.

## demultiplex

Demultiplexation Fastq sequences tool Usage Demultiplexation: %prog -b
[BarCode_file.inp] -i [deep_seq_file.fastq] -o [folder_name] -l 54 -m QUICK
--misreads_cutoff_cons 2

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

## trimming

Trimming Fastq sequences tool Usage Trimming: %prog -d [demultiplexed
Folder]-b [BarCode_file.inp] -q [Quality threshold] -m [method] --output_fmt
fasta

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

## oligo library generator
