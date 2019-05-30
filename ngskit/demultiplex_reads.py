
from __future__ import print_function
from __future__ import unicode_literals
import os
import sys
import pandas as pd
import Levenshtein as Leven
import logging
import argparse
import time
import gzip

import  ngskit.barcodes as barcodes
import  ngskit.qc as qc
#import qc
#import barcodes

class Output_agent(object):
    """Managment of the output demultiplexation

        Parameters
        ----------
        barcodes_list : array_like
            array witht he barcode objects

        out_dir : str
            output directory

        Attributes
        ----------
        output_handlers: dict
             dictionary with all the file handlers
        output_dir: str
            output directory


    """
    def __init__(self, barcodes_list, out_dir):

        self.output_handlers = dict()
        self.out_dir = out_dir
        # open barcode file handlers for barcode
        # improve I/O
        for barcode in barcodes_list:
            self.add_output(barcode)


    def add_output(self, barcode, read=False):


        if read:
            file_name = barcode.id + "_" + str(read['target_len']) + "_F.fastq"
            fileid = barcode.id + "_" + str(read['target_len'])
        else:
            file_name = barcode.id + "_" + str(barcode.trgt_len) + "_F.fastq"
            fileid = barcode.id + "_" + str(barcode.trgt_len)

        out_path = os.path.join(self.out_dir, barcode.id, file_name).replace('\\', '/')

        self.output_handlers[fileid] = open(out_path, 'a')


    def save_seq(self, read1_id, read1_seq, read1_qual, barcode, read):
        """Save sequences in a Fastq file. ToDO: use fastq_tools to improve I/O

        Parameters
        ----------

        Returns
        -------

        """
        fileid = barcode.id + "_" + str(read['target_len'])
        if fileid not in  self.output_handlers:
            self.add_output(barcode, read)

        # save seq in fastq format
        for  _data in [read1_id, read1_seq,"+\n" ,read1_qual]:
            if type(_data) == bytes:
                _data = _data.decode("utf-8")

            self.output_handlers[fileid].write(_data)
#        self.output_handlers[fileid].write(read1_seq)
#        self.output_handlers[fileid].write("+\n")
#        self.output_handlers[fileid].write(read1_qual)


    def close(self):
        for file_handler in self.output_handlers.values():
            file_handler.close()

        return


def match(seq, target, cutoff):
    '''Method used to compare sequence, such barcodes and constant regions.

    Parameters
    ----------
    seq : str
        Nucleotide Sequence trimed from the fastq file .
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
    if type(seq) == bytes:
       seq = seq.decode("utf-8")
       
    distance = count_mismatches(seq, target)

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
        Nucleotide Sequence trimed from the fastq file .
    target : str
        Template to compare, barcode or constant region.

    Returns
    -------
    int
        Mismatches
    '''
    distance = Leven.distance(seq, target)

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
        `dynamic` like `standard` but insertions and deletions in target
         sequence are allowed

    Returns
    -------
        function


    Raises
    ------
        Value Error, if the method to do not exist


    """
    logger = logging.getLogger(__name__)
    
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

    def dynamic(read_seq, barcode, over_end=10, **Kargs):
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
                # Control dynamic length, speed up on assambled Reverse&Forward
                # and diff population of sequences
                if over_end + barcode.b1_len + barcode.c1_len > len(read_seq):
                    return read_map
                # FLOATING WINDOW
                # start searching from the very edge of the cons region 1 (empty vectors)
                # and end  a bit further (insertions)
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

        return read_map

    if method == 'quick':
        return quick
    elif method == 'standard':
        return standard
    elif method == 'simple':
        return simple
    elif method == 'dynamic':
        return dynamic
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
        gbl_stats = qc.Stats(barcodes_list)
    # Temporal, to finish
    dump = Kargs.get('dump', False)

    # open barcode file handlers
    output = Output_agent(barcodes_list, out_dir)

    # Q&D hack, need some elegan solucion for this in the future
    #this can easy fail
    if 'fastq.gz' in inputfile:
        o = gzip.open
    else:
        o = open

    # Open Forward FastaQ file
    with o(inputfile, 'r') as read1:

        for read1_id in read1:
            # Read 4 by 4
            # ID lane info, seq info etc
            # Read seq and Quality info

            read1_seq, read1_strand, read1_qual = [next(read1) for _ in range(3)]

            # For each barcode
            for barcode in barcodes_list:

                read_match_info = identify_seq(read1_seq, barcode, **Kargs)

                if read_match_info['map']:
                    #TODO add some cache to reduce i/o
                    output.save_seq(read1_id, read1_seq, read1_qual,
                             barcode, read_match_info)
                if dump:
                    pass
                if save_frequencies:
                    # save
                    gbl_stats.write_freqs(read_match_info, barcode)
                    pass
    # close
    output.close()
    #move this out of here
    if save_frequencies:
        # write file
        time_stamp = time.ctime()
        fastq_file_name = os.path.basename(inputfile)
        time_out = '{4}_{1}_{2}_{0}_{3}'.format(*time_stamp.split())
        fstats_name =  out_dir + '/Stats/Stats_'+ fastq_file_name + '_' + out_dir + '_'+ Kargs['barcode_file'] + time_out
        gbl_stats.save(fstats_name)

    return





def dump(self, option):
    """Dump information.

    Parameters
    ----------

    Returns
    -------

    """
    pass
    return


def makeoutputdirs(barcode_list, output_dir, is_barcode=True):
    """Helper function to create barcode output files if they not exist.

    Parameters
    ----------

    barcode_list : list,iterable
        list of barcode objects
    output_dir : str
        Path of the output data, by default demultiplex


    """
    logger = logging.getLogger(__name__)

    for sample in barcode_list:
        # Make folders for each barcode
        if is_barcode:
            folder_name = sample.id 
        else:
            folder_name = sample

        out_folder = os.path.join(output_dir, folder_name).replace('\\', '/')
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

    parser.add_argument('-i', '--input_fastqs', nargs='+',
                        dest="input_fastqs", default=False, help='input_fastqs \
                        FASTQ file or files (demultiplex)', required=True)

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
                                 'dynamic'],
                        help="""Type of demultiplexation by default; STANDARD \n
                        `quick`: Only the first barcode and constant region
                        will be  check \n
                        `standard`: Both barcodes and constant regions will be
                         check\n
                        `simple`: Only the barcodes are used \n
                        `dynamic`: frame shift search, Flexible search of\
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


def main():
    """Pipeline Control.

    Parameters
    ----------
    opts


    """

    # Read argtments
    opts = get_options()
    folders_list =  ['Stats', 'Logs']
    makeoutputdirs(folders_list, opts.out_dir, is_barcode=False)

    # init logging
    time_stamp = time.ctime()
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename= opts.out_dir + '/Logs/' +'Dmultplx_'+opts.out_dir+'_'+opts.barcode_file+'_{4}_{1}_{2}_{0}_{3}.log'.format(*time_stamp.split()),
                        filemode='w')

    logger = logging.getLogger(__name__)

    logger.info('JOB START {4} {1} {2} {0} {3}'.format(*time_stamp.split()))


    # Check inputs
    # FASTAQ
    fastqs = opts.input_fastqs
    for fastq in fastqs:
        try:
            os.path.isfile(fastq)
        except:
            raise ValueError('FastQ input {} do not exist'.format(fastq))

    unique = set(fastqs)
    if len(unique) != len(fastqs):
        raise ValueError('duplicate input files')
    # Load Barcodes info
    # check barcodes integrity, peplength, fastq
    barcodes_list = barcodes.read(opts.barcode_file)
    # make output folder
    folders_list = barcodes_list + ['Stats', 'Logs']
    makeoutputdirs(barcodes_list, opts.out_dir)

    # Init Logging
    logger.info('#### DEMULTIPLEXING ####')
    logger.info('Method: {}'.format(opts.dpx_method))
    logger.info('FastQ: {}'.format(opts.input_fastqs))
    logger.info('Barcode: {}'.format(opts.barcode_file))
    # logger.info('Target: {}'.format(opts.target_len))

    logger.info('Misreadings_Barcode: {}'.format(opts.misreads_cutoff_barcode))
    logger.info('Misreadings_Constant: {}'.format(opts.misreads_cutoff_cons))
    logger.info('Stats: {}'.format(opts.save_frequencies))

    # Call to the action
    # To do,: Allow pair end reads.
    for fastq in fastqs:
        logger.info('working on: %s',fastq )
        single_end(fastq, barcodes_list, **opts.__dict__)
        logger.info('next...' )



    # DONE
    time_stamp = time.ctime()
    logger.info('JOB ENDS {4} {1} {2} {0} {3}'.format(*time_stamp.split()))


    return


if __name__ == '__main__':
    main()
