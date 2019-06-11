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


class Output_manager(object):
    """Managment of the output data during demultiplexation

        Attributes
        ----------
        output_handlers: dict
             dictionary with all the file handlers
        output_dir: str
            output directory


    """


    def __init__(self, barcodes_list, out_dir):
        """Init method of the class

        Parameters
        ----------
        barcodes_list : array_like
            array witht he barcode objects

        out_dir : str
            output directory
        """
        self.output_handlers = dict()
        self.out_dir = out_dir

        # open barcode file handlers for barcode
        # improve I/O
        for barcode in barcodes_list:
            self.add_output(barcode)


    def add_output(self, barcode, read=False):
        """Add output in their sample specific file

        Parameters
        ----------
        barcodes_list : array_like
            array witht he barcode objects

        read : bool
            output directory
        """

        if read:
            file_name = barcode.id + "_" + str(read['target_len']) + "_F.fastq"
            fileid = barcode.id + "_" + str(read['target_len'])
        else:
            file_name = barcode.id + "_" + str(barcode.trgt_len) + "_F.fastq"
            fileid = barcode.id + "_" + str(barcode.trgt_len)

        out_path = os.path.join(self.out_dir, barcode.id, file_name).replace('\\', '/')

        self.output_handlers[fileid] = open(out_path, 'a')


    def save_seq(self, read1_id, read1_seq, read1_qual, barcode, read):
        """Save sequences in a Fastq file. 

        Parameters
        ----------
            read1_id: str
                Head of the sequence
            read1_seq: str
                nucleotide seq
            read1_qual: str
                Phred quality string
            barcode: object
                Barcode information
            read; dict
                details of the match

        """
        fileid = barcode.id + "_" + str(read['target_len'])
        if fileid not in  self.output_handlers:
            self.add_output(barcode, read)

        # save seq in fastq format
        for  _data in [read1_id, read1_seq,"+\n" ,read1_qual]:
            if type(_data) == bytes:
                _data = _data.decode("utf-8")

            self.output_handlers[fileid].write(_data)


        return

    def close(self):
        """Close Outputfiles, emty buffers
        
        """ 

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



# TO DO, transfor this to a class
class Demultiplexation_method(object):
    """Sequence identification method.

    Attributes
    ----------
    Function to identify sequences in the fastq poll
    standard: (default) sequence match both barcodes and constant regions
    quick: sequence match barcode and constant  region 1
    simple: sequences with both barcodes
    dynamic: like `standard` but insertions and deletions in target
         sequence are allowed
    misreads_cutoff_barcode :  int
                number of mismatches allowed in the barcode regions (default 1)
    misreads_cutoff_cons : int
                number of mismatches allowed in the constant regions (default 1)


    """
    def __init__(self, options):
        """Init the Class.

            Parameters
            ----------
            options: dict
                Dictionary containing configuration paramaters:

            misreads_cutoff_barcode :  int
                number of mismatches allowed in the barcode regions (default 1)
            misreads_cutoff_cons : int
                number of mismatches allowed in the constant regions (default 1)


        """
        logger = logging.getLogger(__name__)
        
        self.cutoff_barcode = options.get('misreads_cutoff_barcode', 1)
        self.cutoff_cons = options.get('misreads_cutoff_cons', 1)
        self.further_end = options.get('further_end', 10)
        self.stride = options.get('stride', 10)
        

 
    
    def _init_results(self, target_len):

        return {'b1': '',
                'c1': '',
                'c2': '',
                'b2': '',
                'map': False,
                'target_len': target_len}


    def quick(self, read_seq, barcode):
        """Sequence selected if Barcode and constant region 1 are there.

        Parameters
        ----------
        read_seq : str
            sequence
        barcode : object
            barcode object with demultiplexation information

        Returns
        -------
        dict
            dictionary with identification information

        """
        # init return dict
        read_map = self._init_results(barcode.trgt_len)



        # Map Bacode 1, always start at pos 0
        b1 = read_seq[0:barcode.b1_len]

        # Check barcode 1
        if match(b1, barcode.b1_seq, self.cutoff_barcode):
            read_map['b1'] = b1
            # Map  constant region 1, next to barcode 1
            cons1 = read_seq[barcode.b1_len:barcode.c1_len + barcode.b1_len]
            # Check constant region 1
            if match(cons1, barcode.c1_seq, self.cutoff_cons):
                # Success return
                read_map['c1'] = cons1
                read_map['map'] = True

        return read_map

    def standard(self, read_seq, barcode):
        """Selected if  Both Forward and reverse Barcodes and constant regions
        are in the seq, without deletions or insertions.

        Parameters
        ----------
        read_seq : str
            sequence
        barcode : object
            barcode object with demultiplexation information

        Returns
        -------
        dict
            dictionary with identification information

        """
        # init return dict
        read_map = self._init_results(barcode.trgt_len)



        # Map Bacode 1, always start at pos 0
        b1 = read_seq[0:barcode.b1_len]

        # Check barcode 1
        if match(b1, barcode.b1_seq, self.cutoff_barcode):
            read_map['b1'] = b1
            # Extract constant region 1, next to barcode 1
            cons1 = read_seq[barcode.b1_len:barcode.c1_len + barcode.b1_len]
            # Check constantant region 1
            if match(cons1, barcode.c1_seq, self.cutoff_cons):
                # map constant region 2
                read_map['c1'] = cons1
                cons2 = read_seq[barcode.b1_len + barcode.c1_len +
                                 barcode.trgt_len:
                                 barcode.b1_len + barcode.c1_len +
                                 barcode.trgt_len + barcode.c2_len]
                # Check constantant region 1
                if match(cons2, barcode.c2_seq, self.cutoff_cons):
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
                    if match(b2, barcode.b2_seq, self.cutoff_barcode):
                        # Success return

                        read_map['b2'] = b2
                        read_map['map'] = True


        return read_map

    #TODO rename this method
    def simple(self, read_seq, barcode):
        """if Forward and reverse Barcodes are in the seq.
        without deletions or insertions.

        Parameters
        ----------
        read_seq : str
            sequence
        barcode : object
            barcode object with demultiplexation information

        Returns
        -------
        dict
            dictionary with identification information

        """
        # Init results dict
        # init return dict
        read_map = self._init_results(barcode.trgt_len)

   

        # Map Bacode 1, always start at pos 0
        b1 = read_seq[0:barcode.b1_len]
        # Check barcode  1
        if match(b1, barcode.b1_seq, self.cutoff_barcode):
            read_map['b1'] = b1
            # Map Barcode 2
            b2 = read_seq[barcode.b1_len + barcode.c1_len + barcode.trgt_len +
                          barcode.c2_len:
                          barcode.b1_len + barcode.c1_len + barcode.trgt_len +
                          barcode.c2_len + barcode.b2_len]
            # Check Barcode 2
            if match(b2, barcode.b2_seq, self.cutoff_barcode):
                # success
                read_map['b2'] = b2
                read_map['map'] = True


        return read_map

    def dynamic(self, read_seq, barcode):
        """if Forward and reverse Barcodes and constant regions  are in the seq.
        deletions or insertions in the target sequences are allowed.

        Parameters
        ----------
        read_seq : str
            sequence
        barcode : object
            barcode object with demultiplexation information
        further_end : int
            How many up stream over the theorical end of the target position
            should I look for the second constant region

        Returns
        -------
        dict
            dictionary with identification information


        """
        # init return dict
        read_map = self._init_results(barcode.trgt_len)

        # Map Bacode 1, always start at pos 0
        b1 = read_seq[0:barcode.b1_len]

        # Check barcode 1
        if match(b1, barcode.b1_seq, self.cutoff_barcode):
            read_map['b1'] = b1
            # Map constant region 1, next to barcode 1
            cons1 = read_seq[barcode.b1_len:barcode.c1_len + barcode.b1_len]
            # Check constant region
            if match(cons1, barcode.c1_seq, self.cutoff_cons):
                read_map['c1'] = cons1

                # add target distance to self.futher_end, how far up stream
                # it will check
                # (insertions)
                self.further_end += barcode.trgt_len
                # Control dynamic length, speed up on assambled Reverse&Forward
                # and diff population of sequences
                if self.further_end + barcode.b1_len + barcode.c1_len > len(read_seq):
                    return read_map
                # FLOATING WINDOW
                # start searching from the very edge of the cons region 1 (empty vectors)
                # and end  a bit further (insertions)
                for var_target_len in range(self.further_end):

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
                             self.cutoff_cons):
                        # save results
                        read_map['c2'] = dynamic_cons2
                        read_map['target_len'] = var_target_len

                        b2 = read_seq[barcode.b1_len + barcode.c1_len +
                                      var_target_len + barcode.c2_len:
                                      barcode.b1_len + barcode.c1_len +
                                      var_target_len + barcode.c2_len +
                                      barcode.b2_len]
                        # Check Barcode 2
                        if match(b2, barcode.b2_seq, self.cutoff_barcode):
                            read_map['b2'] = b2
                            read_map['map'] = True
                            return read_map

        return read_map


def single_end(inputfile, barcodes_list, out_dir, dpx_method, options):
    """.Main method to demultiplex Single end sequences from deepSeq.
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
    logger = logging.getLogger(__name__)

    # init class
    demultiplexation = Demultiplexation_method(options)

    # Select  method    
    try:
        identify_seq = getattr(demultiplexation, dpx_method.lower(),)
    except AttributeError:
        logger.Warning('Methods {} does not exit using STANDARD'.format(dpx_method))
        identify_seq =  getattr(demultiplexation, 'standard',)
    
    # Control variables
    save_frequencies = options.get('save_frequencies', True)

    if save_frequencies:
        # create pandas file to save stats or go to the log file
        gbl_stats = qc.Stats(barcodes_list)
    # Temporal, to finish
    dump = options.get('dump', False)

    # open barcode file handlers
    output = Output_manager(barcodes_list, out_dir)

    # Q&D hack, need some elegan solucion for this in the future
    #this can easy fail
    if 'fastq.gz' in inputfile:
        o = gzip.open
    else:
        o = open

    # Open Forward FastaQ file
    # eventually I should drop py27 support. 
    with o(inputfile, 'rt', encoding='utf-8') as read1:

        for read1_id in read1:
            # Read 4 by 4
            # ID lane info, seq info etc
            # Read seq and Quality info

            read1_seq, read1_strand, read1_qual = [next(read1) for _ in range(3)]

            # For each barcode
            for barcode in barcodes_list:

                read_match_info = identify_seq(read1_seq, barcode)

                if read_match_info['map']:
                                        
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
        seconds_time = int(time.time())
        fastq_filename = os.path.basename(inputfile)
        barcode_filename = options['barcode_file']
        fstats_name =  f'{out_dir}/Stats/Stats_{fastq_filename}_{out_dir}_{barcode_filename}_{seconds_time}'
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


def makeoutputdirs(barcode_list, output_dir):
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
        if type(sample)!= str:
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

def check_fastq_exist(opts):
    """check all fastq Files exist.

    Parameters
    ----------
        options
    Returns
    -------

    """

    fastqs = opts.input_fastqs
    for fastq in fastqs:
        try:
            os.path.isfile(fastq)
        except:
            raise ValueError('FastQ input {} do not exist'.format(fastq))

    return


def check_fastq_unique(opts):
    # Files are unique
    
    fastqs = opts.input_fastqs
    unique = set(fastqs)
    if len(unique) != len(fastqs):
        raise ValueError('duplicate input files')

    return

def checking_inputs(opts):
    """wrapper of control of input files

    Parameters
    ----------
        options

    Returns
    -------
        error if a file do not complain

    """

    check_fastq_exist(opts)
    check_fastq_unique(opts)
    
    # move this to barcodes
    # Load Barcodes info
    # check barcodes integrity, peplength, fastq
    barcodes_list = barcodes.read(opts.barcode_file)
    [bc.sanity_check() for bc in barcodes_list]

    return 


##############################################################
##############################################################



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

    parser.add_argument('--further_end', action="store",
                        dest="further_end", default=10, type=int,
                        help='Number of positions to keep looking \
                        for the 2nd constant region after the target lenght. \
                        Only applys on dynamic method')

    parser.add_argument('--dump', help='Dump constant regions', dest='dump',
                        default=False, action='store_true')

    parser.add_argument('--no-save_frequencies', help='Do not Save match \
                        frequencies', dest='save_frequencies', default=True,
                        action='store_false')

    options = parser.parse_args()

    return options


def main():
    """Main Pipeline.

    Parameters
    ----------
    opts


    """

    # Read arguments
    opts = get_options()

    # init logging
    seconds_time = int(time.time())
    time_stamp = time.ctime()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=  f'{opts.out_dir}/Logs/Dmultplx_{opts.out_dir}_{opts.barcode_file}_{seconds_time}.log',
                        filemode='w')

    logger = logging.getLogger(__name__)

    logger.info('JOB START {4} {1} {2} {0} {3}'.format(*time_stamp.split()))

    # Check inputs consistency
    checking_inputs(opts)

    # init barcodes
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
    fastqs = opts.input_fastqs
    for fastq in fastqs:
        logger.info('working on: %s',fastq )
        single_end(fastq, barcodes_list, opts.out_dir, opts.dpx_method, opts.__dict__)
        logger.info('next...' )



    # DONE
    time_stamp = time.ctime()
    logger.info('JOB ENDS {4} {1} {2} {0} {3}'.format(*time_stamp.split()))


    return


if __name__ == '__main__':
    main()
