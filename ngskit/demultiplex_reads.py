from __future__ import print_function
import os
import sys
import pandas as pd
import Levenshtein as leven
import logging
import argparse
import time

import barcodes


def match(seq, target, cutoff):
    '''Method used to compare sequence, such barcodes and constant regions.

    Parameters
    ----------
    seq : str
        Nucleotide Sequence trimed from the fastaq file .
    target : str
        Template to compare, barcode or constant region.
    cutoff : str
        Maximum number of diff accepted between the two sequnces.
    read1_id : str
        Id of the sequences, to debug purposes. This info and the number of
        differences is recorder in the log file

    Returns
    -------
    boolean
        False if there are diff over the cutoff

    '''
    cutoff = int(cutoff)
    distance = leven.distance(seq, target)

    if distance > cutoff:

        return False
    else:

        return True


def count_mismatches(seq, target):
    '''Method used to compare sequence.
    such barcodes and constant regions.

    Parameters
    ----------
    seq : str
        Nucleotide Sequence trimed from the fastaq file .
    target : str
        Template to compare, barcode or constant region.

    Returns
    -------
    int
        Mismatches
    '''
    distance = leven.distance(seq, target)

    return distance


def identification_method(method='standard'):
    """Sequence identification method.

    Parameters
    ----------
    method : str
        Fucntion to identify ssequences in the fastq poll
        `standard` (default) sequence match both barcodes and constant regions
        `quick` sequence match barcode and constant  region 1
        `simple` sequences with both barcodes
        `float_window` like `standard` but insertions and deletions in target
         sequence are allowed

    Returns
    -------
        function


    Raises
    ------
        Value Error, if the method to do not exist


    """

    def quick(read_seq, barcode, **Kargs):
        """Sequence selected if Barcode and constant region 1 are there.

        Parameters
        ----------
        read_seq : str
            sequence
        barcode : object
            barcode object with demultiplexation information
        misreads_cutoff_barcode :  int
            number of mismatches allowed in the barcode regions (default 1)
        misreads_cutoff_cons : int
            number of mismatches allowed in the constant regions (default 1)
        Returns
        -------
        dict
            dictionary with identification information

        """
        # init return dict
        read_map = {'b1': '',
                    'c1': '',
                    'c2': '',
                    'b2': '',
                    'map': False,
                    'target_len': barcode.trgt_len}
        # get mismatch configuration
        # try:
        misreads_cutoff_barcode = Kargs['misreads_cutoff_barcode']
        misreads_cutoff_cons = Kargs['misreads_cutoff_cons']
        # except KeyError:
        #     # logger.info('Misreads configuration missing using defaults')
        # misreads_cutoff_barcode = Kargs.get('misreads_cutoff_barcode', 0)
        # misreads_cutoff_cons = Kargs.get('misreads_cutoff_cons', 2)

        # Map Bacode 1, always start at pos 0
        b1 = read_seq[0:barcode.b1_len]

        # Check barcode 1
        if match(b1, barcode.b1_seq, misreads_cutoff_barcode):
            read_map['b1'] = b1
            # Map  constant region 1, next to barcode 1
            cons1 = read_seq[barcode.b1_len:barcode.c1_len + barcode.b1_len]
            # Check constant region 1
            if match(cons1, barcode.c1_seq, misreads_cutoff_cons):
                # Success return
                read_map['c1'] = cons1
                read_map['map'] = True

        return read_map

    def standard(read_seq, barcode, **Kargs):
        """Selected if  Both Forward and reverse Barcodes and constant regions
        are in the seq, without deletions or insertions.

        Parameters
        ----------
        read_seq : str
            sequence
        barcode : object
            barcode object with demultiplexation information
        misreads_cutoff_barcode :  int
            number of mismatches allowed in the barcode regions (default 1)
        misreads_cutoff_cons : int
            number of mismatches allowed in the constant regions (default 1)

        Returns
        -------
        dict
            dictionary with identification information

        """
        # Init results dict
        read_map = {'b1': '',
                    'c1': '',
                    'c2': '',
                    'b2': '',
                    'map': False,
                    'target_len': barcode.trgt_len}

        misreads_cutoff_barcode = Kargs['misreads_cutoff_barcode']
        misreads_cutoff_cons = Kargs['misreads_cutoff_cons']

        # Map Bacode 1, always start at pos 0
        b1 = read_seq[0:barcode.b1_len]

        # Check barcode 1
        if match(b1, barcode.b1_seq, misreads_cutoff_barcode):
            read_map['b1'] = b1
            # Extract constant region 1, next to barcode 1
            cons1 = read_seq[barcode.b1_len:barcode.c1_len + barcode.b1_len]
            # Check constantant region 1
            if match(cons1, barcode.c1_seq, misreads_cutoff_cons):
                # map constant region 2
                read_map['c1'] = cons1
                cons2 = read_seq[barcode.b1_len + barcode.c1_len +
                                 barcode.trgt_len:
                                 barcode.b1_len + barcode.c1_len +
                                 barcode.trgt_len + barcode.c2_len]
                # Check constantant region 1
                if match(cons2, barcode.c2_seq, misreads_cutoff_cons):
                    read_map['c2'] = cons2
                    # Map design, target seq
                    # design = read_seq[barcode.b1_len + barcode.c1_len:
                    #                   barcode.b1_len + barcode.c1_len +
                    #                   barcode.trgt_len]
                    # Map barcode 2
                    b2 = read_seq[barcode.b1_len + barcode.c1_len +
                                  barcode.trgt_len + barcode.c2_len:
                                  barcode.b1_len + barcode.c1_len +
                                  barcode.trgt_len + barcode.c2_len +
                                  barcode.b2_len]
                    # Check Barcode 2
                    if match(b2, barcode.b2_seq, misreads_cutoff_barcode):
                        # Success return

                        read_map['b2'] = b2
                        read_map['map'] = True

        return read_map

    def simple(read_seq, barcode, **Kargs):
        """if Forward and reverse Barcodes are in the seq.
        without deletions or insertions.

        Parameters
        ----------
        read_seq : str
            sequence
        barcode : object
            barcode object with demultiplexation information
        misreads_cutoff_barcode :  int
            number of mismatches allowed in the barcode regions (default 1)

        Returns
        -------
        dict
            dictionary with identification information

        """
        # Init results dict
        read_map = {'b1': '',
                    'c1': '',
                    'c2': '',
                    'b2': '',
                    'map': False,
                    'target_len': barcode.trgt_len}

        misreads_cutoff_barcode = Kargs['misreads_cutoff_barcode']

        # Map Bacode 1, always start at pos 0
        b1 = read_seq[0:barcode.b1_len]
        # Check barcode  1
        if match(b1, barcode.b1_seq, misreads_cutoff_barcode):
            read_map['b1'] = b1
            # Map Barcode 2
            b2 = read_seq[barcode.b1_len + barcode.c1_len + barcode.trgt_len +
                          barcode.c2_len:
                          barcode.b1_len + barcode.c1_len + barcode.trgt_len +
                          barcode.c2_len + barcode.b2_len]
            # Check Barcode 2
            if match(b2, barcode.b2_seq, misreads_cutoff_barcode):
                # success
                read_map['b2'] = b2
                read_map['map'] = True

        return read_map

    def float_window(read_seq, barcode, over_end=10, **Kargs):
        """if Forward and reverse Barcodes and constant regions  are in the seq.
        deletions or insertions in the target sequences are allowed.

        Parameters
        ----------
        read_seq : str
            sequence
        barcode : object
            barcode object with demultiplexation information
        misreads_cutoff_barcode :  int
            number of mismatches allowed in the barcode regions (default 1)
        misreads_cutoff_cons : int
            number of mismatches allowed in the constant regions (default 1)
        over_end : int
            How many up stream over the theorical end of the target position
            should I look for the second constant region

        Returns
        -------
        dict
            dictionary with identification information


        """
        # Init results dict
        read_map = {'b1': '',
                    'c1': '',
                    'c2': '',
                    'b2': '',
                    'map': False,
                    'target_len': barcode.trgt_len}

        misreads_cutoff_barcode = Kargs['misreads_cutoff_barcode']
        misreads_cutoff_cons = Kargs['misreads_cutoff_cons']

        # Map Bacode 1, always start at pos 0
        b1 = read_seq[0:barcode.b1_len]

        # Check barcode 1
        if match(b1, barcode.b1_seq, misreads_cutoff_barcode):
            read_map['b1'] = b1
            # Map constant region 1, next to barcode 1
            cons1 = read_seq[barcode.b1_len:barcode.c1_len + barcode.b1_len]
            # Check constant region
            if match(cons1, barcode.c1_seq, misreads_cutoff_cons):
                read_map['c1'] = cons1

                # add target distance to over_end, how far up stream
                # it will check
                # (insertions)
                over_end += barcode.trgt_len

                # FLOATING WINDOW
                # start from the very edge of the cons region 1 (empty vectors)
                for var_target_len in range(over_end):

                    # extract constant region 2, using float window
                    dynamic_cons2 = read_seq[barcode.b1_len + barcode.c1_len +
                                             var_target_len: barcode.b1_len +
                                             barcode.c1_len + var_target_len +
                                             barcode.c2_len]
                    logger.debug('{}_{}_{}:{}_{}_{}_{}'.format(barcode.b1_len,
                                                               barcode.c1_len,
                                                               var_target_len,
                                                               barcode.b1_len,
                                                               barcode.c1_len,
                                                               var_target_len,
                                                               barcode.c2_len))
                    logger.debug(dynamic_cons2)
                    # check dynamic region against constant region 2
                    if match(dynamic_cons2, barcode.c2_seq,
                             misreads_cutoff_cons):
                        # save results
                        read_map['c2'] = dynamic_cons2
                        read_map['target_len'] = var_target_len

                        b2 = read_seq[barcode.b1_len + barcode.c1_len +
                                      var_target_len + barcode.c2_len:
                                      barcode.b1_len + barcode.c1_len +
                                      var_target_len + barcode.c2_len +
                                      barcode.b2_len]
                        # Check Barcode 2
                        if match(b2, barcode.b2_seq, misreads_cutoff_barcode):
                            read_map['b2'] = b2
                            read_map['map'] = True

        return read_map

    if method == 'quick':
        return quick
    elif method == 'standard':
        return standard
    elif method == 'simple':
        return simple
    elif method == 'float_window':
        return float_window
    else:
        raise ValueError('Method {} do not exist'.format(method))


def single_end(inputfile, barcodes_list, out_dir='demultiplex',
               dpx_method='standard', **Kargs):
    """.

    Main method to demultiplex Single end sequences from deepSeq.
    No quality score applyed [need improvements]


    Parameters
    ----------
        inputfile : str
            Fastq file
        barcodes : list
            list or iterable with barcode objects
        out_dir : str
            path to save demultiplexed fastq files

        dpx_type : str
            Type of demultiplexation by default `standard`


    See Also
    --------
    identification_method

    """
    # Select  method
    identify_seq = identification_method(dpx_method)
    # Control variables
    save_frequencies = Kargs.get('save_frequencies', True)

    if save_frequencies:
        # create pandas file to save stats or go to the log file
        gbl_stats = init_freqs_track(barcodes_list)
    # Temporal, to finish
    dump = Kargs.get('dump', False)

    # Open Forward FastaQ file
    with open(inputfile, 'r') as read1:

        for read1_id in read1:
            # Read 4 by 4
            # ID lane info, seq info etc
            # Read seq and Quality info

            read1_seq, read1_strand, read1_qual = [next(read1) for _ in range(3)]

            # For each barcode
            for barcode in barcodes_list:

                read_match_info = identify_seq(read1_seq, barcode, **Kargs)

                if read_match_info['map']:
                    save_seq(read1_id, read1_seq, read1_qual,
                             barcode, read_match_info, out_dir)
                if dump:
                    pass
                if save_frequencies:
                    # save
                    gbl_stats = write_freqs(read_match_info, barcode, gbl_stats)
                    pass
    # close
    if save_frequencies:
        # write file
        time_stamp = time.ctime()
        gbl_stats.to_csv(out_dir + '_'+ Kargs['barcode_file'] +'_stats_{4}_{1}_{2}_{0}_{3}.csv'.format(*time_stamp.split()))

    return


def save_seq(read1_id, read1_seq, read1_qual, barcode, read, out_dir):
    """Save sequences in a Fastq file. ToDO: use fastq_tools to improve I/O

    Parameters
    ----------

    Returns
    -------

    """
    file_name = barcode.id + "_" + str(read['target_len']) + "_F.fastq"
    out_file = os.path.join(out_dir, barcode.id, file_name).replace('\\', '/')
    f = open(out_file, 'a')
    f.write(read1_id)
    f.write(read1_seq)
    f.write("+\n")
    f.write(read1_qual)
    f.close()

    return


def init_freqs_track(barcodes_list):
    """Init dataframe to Save frequencies.

    Parameters
    ----------
    barcodes : list, iterable
        contains barcode objects to demultiplex sequneces

    Returns
    -------
    DataFrame

    idx       b1  c1  c2  b2 target_len1 ... target_lenX
    sample1
    sample2
    ...
    sampleX

    """
    index_ids = []
    columns_id = ['b1', 'c1', 'c2', 'b2']
    for barcode in barcodes_list:

        index_ids.append(barcode.id)
        if not str(barcode.trgt_len) in columns_id:
            columns_id.append(str(barcode.trgt_len))

    # init stats

    gbl_stats = pd.DataFrame(index=index_ids, columns=columns_id)
    gbl_stats.fillna(0, inplace=True)
    logger.debug(gbl_stats)

    return gbl_stats


def write_freqs(read, barcode, gbl_stats):
    """Update frequencies in the DataFrame.

    Parameters
    ----------
    read : dict
        Dict from the identification fuction with information about the match
    barcode : object

    Returns
    -------
    DataFrame

    """
    if read['b1'] != '':
        # Seq match barcode 1 (+1)
        gbl_stats.at[barcode.id, 'b1'] = gbl_stats.at[barcode.id,
                                                      'b1'] + 1
    if read['c1'] != '':
        # Seq match conc 1 (+1)
        gbl_stats.at[barcode.id, 'c1'] = gbl_stats.at[barcode.id,
                                                      'c1'] + 1
    if read['c2'] != '':
        # Seq match conc 2 (+1)
        gbl_stats.at[barcode.id, 'c2'] = gbl_stats.at[barcode.id,
                                                      'c2'] + 1
    if read['b2'] != '':
        # Seq match barcode 2 (+1)
        gbl_stats.at[barcode.id, 'b2'] = gbl_stats.at[barcode.id,
                                                      'b2'] + 1
    # Target sequnces, length +1
    if read['map']:
        if str(read['target_len']) in gbl_stats.columns:
            gbl_stats.at[barcode.id,
                         str(read['target_len'])] = gbl_stats.at[barcode.id,
                                                                 str(read['target_len'])] + 1

        else:
            gbl_stats.at[barcode.id, str(read['target_len'])] = 1
            gbl_stats.fillna(0, inplace=True)

    return gbl_stats


def dump(self, option):
    """Dump information.

    Parameters
    ----------

    Returns
    -------

    """
    pass
    return


def makeoutputdirs(barcode_list, output_dir):
    """Helper function to create barcode output files if they not exist.

    Parameters
    ----------

    barcode_list : list,iterable
        list of barcode objects
    output_dir : str
        Path of the output data, by default demultiplex


    """
    for sample in barcode_list:
        # Make folders for each barcode
        out_folder = os.path.join(output_dir, sample.id).replace('\\', '/')
        try:
            os.makedirs(out_folder)
        except OSError:
            e = sys.exc_info()
            # print 'Warning, Folder',out_folder,' already exist'
            print('Warning {}'.format(e[1]))
            logger.info('Folder {} already exist'.format(out_folder))

    return

##############################################################
##############################################################
##############################################################
##############################################################
#####################


def get_options():
    """Get arguments from command line.

    Parameters
    ----------

    Returns
    -------

    """
    parser = argparse.ArgumentParser(description="""
    Demultiplexation Fastq sequences tool

    Usage Demultiplexation:
    %prog -b [BarCode_file.inp]  -i [deep_seq_file.fastq] -o [folder_name]\
    -l 54 -m QUICK --misreads_cutoff_cons 2

    """)

    # parser.add_argument('-a', '--action',
    #                     type='choice',
    #                     action='store',
    #                     dest='action',
    #                     choices=['singleEND', 'pairEND'],
    #                     default=False,
    #                     help='to develop')

    # Change this name, get thisn info from the barcodes
    # Depricated, information added in the barcodes file,
    # allow multiple lenghts
    # parser.add_argument('-l', '--target_len', action="store",
    #                     dest="target_len",
    #                     default=False, help='Lengh of the designed oligo',
    #                     required=True)

    parser.add_argument('-i', '--input_fastq', action="store",
                        dest="input_fastq", default=False, help='input_fastq \
                        FASTAQ (demultiplex)', required=True)

    parser.add_argument('-b', '--barcode_file', action="store",
                        dest="barcode_file", default=False, help='File that \
                        contains barcodes and cosntant regions', required=True)

    parser.add_argument('-o', '--out_dir', action="store", dest="out_dir",
                        default='demultiplex', help='Output folder, called \
                        demultiplex by default')

    # optional Arguments
    parser.add_argument('-m', '--demultiplexation_method', action="store",
                        dest="dpx_method", default='standard', type=str,
                        choices=['quick',
                                 'standard',
                                 'simple',
                                 'float_window'],
                        help="""Type of demultiplexation by default; STANDARD \n
                        `quick`: Only the first barcode and constant region
                        will be  check \n
                        `standard`: Both barcodes and constant regions will be
                         check\n
                        `simple`: Only the barcodes are used \n
                        `float_window`: frame shift search, Flexible search of\
                        the second  constant region and barcode\n
                        """)
    # Default 1
    parser.add_argument('--misreads_cutoff_cons', action="store",
                        dest="misreads_cutoff_cons", default=2, type=int,
                        help='Max number of misreading allowed in the constant \
                        constant_region (default 2)')

    parser.add_argument('--misreads_cutoff_barcode', action="store",
                        dest="misreads_cutoff_barcode", default=1, type=int,
                        help='Max number of misreading allowed in the constant \
                        constant_region  (default 1)')

    parser.add_argument('--dump', help='Dump constant regions', dest='dump',
                        default=False, action='store_true')

    parser.add_argument('--no-save_frequencies', help='Do not Save match \
                        frequencies', dest='save_frequencies', default=True,
                        action='store_false')

    options = parser.parse_args()

    return options


def workflow(opts):
    """Pipeline Control.

    Parameters
    ----------
    opts


    """
    # Check inputs
    # FASTAQ
    fastaq = opts.input_fastq
    try:
        os.path.isfile(fastaq)
    except:
        raise ValueError('FastQ input {} do not exist'.format(fastaq))

    # Load Barcodes info
    # check barcodes integrity, peplength, fastq
    barcodes_list = barcodes.read(opts.barcode_file)
    # make output folder
    makeoutputdirs(barcodes_list, opts.out_dir)

    # Init Logging
    logger.info('#### DEMULTIPLEXING ####')
    logger.info('Method: {}'.format(opts.dpx_method))
    logger.info('FastQ: {}'.format(opts.input_fastq))
    logger.info('Barcode: {}'.format(opts.barcode_file))
    # logger.info('Target: {}'.format(opts.target_len))

    logger.info('Misreadings_Barcode: {}'.format(opts.misreads_cutoff_barcode))
    logger.info('Misreadings_Constant: {}'.format(opts.misreads_cutoff_cons))
    logger.info('Stats: {}'.format(opts.save_frequencies))

    # Call to the action
    # To do,: Allow pair end reads.
    single_end(fastaq, barcodes_list, **opts.__dict__)

    return


if __name__ == '__main__':
    # Read argtments
    opts = get_options()

    # init logging
    time_stamp = time.ctime()
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename= 'Demultiplex_'+opts.barcode_file+'_{4}_{1}_{2}_{0}_{3}.log'.format(*time_stamp.split()),
                        filemode='w')
    logger = logging.getLogger(__name__)

    logger.info('JOB START {4} {1} {2} {0} {3}'.format(*time_stamp.split()))
    # DEMULTIPLEX
    workflow(opts)
    # DONE
    time_stamp = time.ctime()
    logger.info('JOB ENDS {4} {1} {2} {0} {3}'.format(*time_stamp.split()))
