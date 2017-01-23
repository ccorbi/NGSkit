from __future__ import print_function
from __future__ import division
import os
import sys
import logging
import argparse
import time

import barcodes
from utils import fasta_tools, fastq_tools



def trimming(demultiplexed_fastq, barcode, quality_threshold,
            trgt_len, output_fmt, output_folder):
    """Extract seq from the FASTAQ demultiplexed files. Trim barcodes + Constant

    Parameters
    ----------
    demultiplexed_fastq : str
        Path of the demultiplexed fastq file
    barcode : barcode.object
        Barcode object wiht info about barcode and constant regions
    quality_threshold : int
        reading quality Threshold, any sequence will be trimmed under that level
    trgt_len : int
        length in bases of the target sequences.
    output_fmt :  str
        Output format, by default fasta
    working_folder : str
        Output folder to save files with trimmed sequences

    Returns
    -------
    output format save fasta or fastq


    Notes
    -----

    Result str, in Fasta format
    >FASTAQ_ID+ length + Quality
    ATGATGGTAGTAGTAGAAAGATAGATGATGATGAT
    it will be storage:
    /data_path/Sequences/Sample_id.fasta


    """
    # Init the output format, retunr a function
    if output_fmt == 'fasta':
        save_seq = fasta_tools.write_fasta_sequence
        filehdl_output = open(output_folder+barcode.id+'.fasta','a')
        logger.info('Output file: %s' % (output_folder+barcode.id+'.fasta'))

    if output_fmt == 'fastq':
        save_seq = fastq_tools.write_fastq_sequence
        filehdl_output = open(output_folder+barcode.id+'.fastq','a')
        logger.info('Output file: %s' % (output_folder+barcode.id+'.fastq'))
    # check barcodes integrity, peplength, fastq
    # barcodes_list = barcodes.read(barcode_file)

    # Stats
    nseqs = 0
    ntrimed = 0
    # Open Fastq file
    with open(demultiplexed_fastq, 'r') as read1:
        for read1_id in read1:
            # Read 4 by 4
            # ID lane info, seq info etc
            # Read seq and Quality info
            read1_seq, read1_strand, read1_qual = [next(read1) for _ in range(3)]
            #Translate the Quality to a list of Integers
            qual = [ord(c)-33 for c in read1_qual.rstrip("\n")]

            target_sequence = read1_seq[barcode.b1_len+barcode.c1_len:
                                        barcode.b1_len+barcode.c1_len+trgt_len]
            #remove the quality of the barcode and the constant region
            target_qual = qual[barcode.b1_len+barcode.c1_len:
                                        barcode.b1_len+barcode.c1_len+trgt_len]
            nseqs += 1
            # Control
            avg_quality = sum(target_qual)/float(len(target_qual))
            if len(target_sequence) == trgt_len and avg_quality >= quality_threshold:

                ntrimed += 1
                # save output format
                # attach Qavg to the id
                seq_id = '{}_{}'.format(read1_id.strip(), avg_quality)
                save_seq([seq_id, target_sequence, target_qual],
                         file_output=filehdl_output)
                # save
            else:
                # Stats
                pass

    logger.info('Read %i Sequences' % (nseqs))
    logger.info('Trimmed %i Sequences' % (ntrimed))
    filehdl_output.close()

def get_options():
    """Get arguments from command line.

    Parameters
    ----------

    Returns
    -------

    """
    parser = argparse.ArgumentParser(description="""
    Trimming Fastq sequences tool

    Usage Trimming:
    %prog -d [demultiplexed Folder]-b [BarCode_file.inp]  -q [Quality threshold]\
    -m [method] --output_fmt fasta

    """)

    parser.add_argument('-d', '--input_folder', action="store",
                        dest="input_folder", default=False, help='Folder \
                        contains demultiplexed folders and files', required=True)

    parser.add_argument('-b', '--barcode_file', action="store",
                        dest="barcode_file", default=False, help='File that \
                        contains barcodes and cosntant regions', required=True)

    parser.add_argument('-o', '--out_folder', action="store", dest="out_folder",
                        default='Sequences', help='Output folder, called \
                        Sequences by default')

    # optional Arguments
    parser.add_argument('-m', '--trimming_method', action="store",
                        dest="trimming_method", default='standard', type=str,
                        choices=['standard',
                                 'float_window'],
                        help="""                        """)
    # Default 1
    parser.add_argument('-q', '--quality', action="store",
                        dest="quality", default=30, type=int,
                        help='Quality reading threshold \
                         (default 30)')


    parser.add_argument('--output_fmt', help='Output format, default fasta',
                        dest='output_fmt', default='fasta', action='store')


    options = parser.parse_args()

    return options

def workflow(opts):
    """Pipeline Control.

    Parameters
    ----------
    opts


    """
    # Check inputs
    # Load Barcodes info
    # check barcodes integrity, peplength, fastq
    barcodes_list = barcodes.read(opts.barcode_file)
    # make output folder
    # Init Logging
    logger.info('#### TRIMMING ####')
    logger.info('Method: {}'.format(opts.trimming_method))
    logger.info('Quality threshold: {}'.format(opts.quality))
    logger.info('Output format: {}'.format(opts.output_fmt))
    #
    logger.info('Barcode file: {}'.format(opts.barcode_file))
    logger.info('Input folder: {}'.format(opts.input_folder))
    output_folder = opts.input_folder+'/'+opts.out_folder+'/'
    logger.info('Output folder: {}'.format(output_folder))

    # Create output folder
    try:
        # by default Sequences
        os.makedirs(output_folder)
    except OSError:
        e = sys.exc_info()
        #print 'Warning, Folder',out_folder,' already exist'
        print('Warning  {}'.format(e[1]))

    # foreach sample in barcodes
    for barcode in barcodes_list:
        logger.info('Triming Sample: {}'.format(barcode.id))
        # folder must == sample id in the barcode
        working_folder = './'+opts.input_folder+'/'+barcode.id+'/'
        # get all fastq under the folder
        for demultiplexed_fastq in os.listdir(working_folder):
            # ToDO: only get fastq files
            #ToDo: only those I want (target lenthg)
            # if method is dynamic, get all the files in the folder
            if opts.trimming_method == 'dynamic':
                # To do
                # raw_name = demultiplexed_file.replace('_F.fastq','')
                # read the length from the file
                pass
            else:
                logger.info('Triming file: {}'.format(demultiplexed_fastq))
                # Trim time
                dir_emultiplexed_fastq = working_folder+demultiplexed_fastq
                trimming(dir_emultiplexed_fastq,
                        barcode,
                        quality_threshold= opts.quality,
                        trgt_len= barcode.trgt_len,
                        output_fmt= opts.output_fmt,
                        output_folder=output_folder)

    return


if __name__ == '__main__':
    # Read argtments
    opts = get_options()

    # init logging
    time_stamp = time.ctime()
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename= 'Trimming_'+opts.input_folder+'_'+opts.barcode_file+'_{4}_{1}_{2}_{0}_{3}.log'.format(*time_stamp.split()),
                        filemode='w')
    logger = logging.getLogger(__name__)

    logger.info('JOB START {4} {1} {2} {0} {3}'.format(*time_stamp.split()))
    # DEMULTIPLEX
    workflow(opts)
    # DONE
    time_stamp = time.ctime()
    logger.info('JOB ENDS {4} {1} {2} {0} {3}'.format(*time_stamp.split()))
