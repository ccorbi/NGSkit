import os
import sys
import barcodes
import logging
import time
import pandas as pd
import argparse


def makeoutputdirs(barcode,out_path):

    '''Helper function to create barcode output files if they not exist

    :param str barcode: barcode name
    :param str output_path: Path of the output data, by default demultiplex'''

    for sample_id in barcode:
        #Make folders for each barcode
        out_folder = os.path.join(out_path, sample_id).replace('\\', '/')
        try:
            os.makedirs(out_folder)
        except OSError:
            e = sys.exc_info()
            #print 'Warning, Folder',out_folder,' already exist'
            print 'Warning', e[1]
    # Open Forward
    return

def read_fasta2dict(filename):

    '''Helper function quickly read a fasta file into a dict [need improvement]

    :param str filename: Fasta file path and name

    :return: Dictionary header is key, sequence is the value
    :rtype: dict '''

    fasta_dic = dict()

    with open(filename,'r') as input_file:
        for line in input_file:
            if line[0] == '>':
                fasta_dic[line[1:].strip()] =''
                last_id = line[1:].strip()
            else:
                fasta_dic[last_id] = line.strip()

    return fasta_dic


    def trim_all_(self,barcode_file,data_path,**Kargs):
        '''
        For each sample in the barcodes file, trim seq from the FASTAQ demultiplexed files
        , right now, just remove the barcodes and filter by qualty and length

        :param str barcode_file: file with the barcodes info to be transformed Dictionary with the barcodes
        :param str data_path: demultiplexing folder '''

        barcodes = read_barcodes(barcode_file)

        for sample_id in barcodes:
            # Get all the files in the sample folder
            for demultiplex_file in os.listdir('./'+data_path+'/'+sample_id):

#                output = demultiplex_file.strip()
                raw_name = demultiplex_file.replace('_F.fastq','')
                hit_len = raw_name.split('_')[-1:]
                #print hit_len,sys.argv[3]
                if hit_len[0]==self.peplength:
                    print demultiplex_file
                    path_dpx_file = './'+data_path+'/'+sample_id+'/'+demultiplex_file
                    self.trim_naive(path_dpx_file,barcodes,sample_id,data_path,**Kargs)

        return

def trim_naive(self,demultiplexed_file,barcodes,sample_id,data_path,**Kargs):

    '''Extract seq from the FASTAQ demultiplexed files,
    remove barcodes and constant region
    filter by quality

    Result str, in Fasta format
    >FASTAQ_ID+ length + Quality
    ATGATGGTAGTAGTAGAAAGATAGATGATGATGAT
    it will be storage:
    /data_path/Sequences/Sample_id.fasta

    :param str demultiplexed_file: Path to the Fastq file demultiplexed
    :param dict barcodes: barcodes info dictoniary
    :param str sample_id: sample_id
    :param str data_path: Folder name where is the demultiplexed info.
    '''
    # self.setup_cfg(**param)


    assert demultiplexed_file!=''
    assert len(barcodes)>0
    assert sample_id!=''

    assert len(data_path) > 0

    assert len(self.peplength) > 0

    chopp = Kargs.get('chopp')


    F_cons = barcodes[sample_id][2]
    S_cons = barcodes[sample_id][3]
    barcode1 = barcodes[sample_id][0]
#        barcode2 = barcodes[sample_id][1]

    result = open("./%s/Sequences/%s.fasta" % (data_path,sample_id), "a")
    no_complaint = 0
    complaint = 0
    with open(demultiplexed_file, 'r') as forward:
        for forward_id in forward:

            header_id = forward_id
            header_id = str(forward_id.partition(' ')[0]).strip()
            sequence, line3, quality = [next(forward) for _ in range(3)]

            forward_id = forward_id.partition(' ')

            #Translate the Quality to numbers
            qual = [ord(c)-33 for c in quality.rstrip("\n")]

            if chopp == 'peptide':
                sequence = sequence[len(barcode1)+len(F_cons):int(self.peplength)+len(barcode1)+len(F_cons)]
                #remove the quality of the barcode and the constant region
                qual = qual[len(barcode1)+len(F_cons):int(self.peplength)+len(barcode1)+len(F_cons)]
                lencutoff = int(self.peplength)

            elif chopp == 'constant':
                sequence = sequence[len(barcode1):int(self.peplength)+len(barcode1)+len(S_cons)+len(F_cons)]
                #remove the quality of the barcode and the constant region
                qual = qual[len(barcode1):int(self.peplength)+len(barcode1)+len(F_cons)+len(S_cons)]
                lencutoff = int(self.peplength)+len(F_cons)+len(S_cons)
            else:
                print 'Chopp method not defined'
                raise KeyError




            #split the remain seq in diff frames in order to find the reverse constant region
            #Not always the seq has the 48bp, that why we use this method.

            dseq_q = sum(qual)/float(len(qual))
            if len(sequence) == lencutoff and dseq_q >= self.reading_quality_cutoff:

                print >> result, ">"+header_id+"_F_"+str(len(sequence))+"_"+str(dseq_q)
                print >> result, sequence
                complaint +=1
            else:
                logger.info('{} {} {} '.format(header_id, dseq_q, sequence))
                no_complaint +=1

    logger.info('{} have been discarted of {} '.format(no_complaint,no_complaint+complaint))

    return



def get_options():
    """Get arguments from command line.

    Parameters
    ----------

    Returns
    -------

    """
    parser = argparse.ArgumentParser(description="""
    Trimming Fastq sequences tool

    Usage Demultiplexation:
    %prog -b [BarCode_file.inp]  -i [deep_seq_file.fastq] -o [folder_name]\
    -l 54 -m QUICK --misreads_cutoff_cons 2

    """)

    parser.add_argument('-i', '--input_fastq', action="store",
                        dest="input_fastq", default=False, help='input_fastq \
                        FASTAQ (demultiplex)')

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
                        dest="misreads_cutoff_cons", default=1, type=int,
                        help='Max number of misreading allowed in the constant \
                        constant_region (default 1)')

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
    # Load Barcodes info
    # check barcodes integrity, peplength, fastq
    barcodes_list = barcodes.read(opts.barcode_file)


    # Init Logging
    logger.info('#### TRIMMING ####')
    logger.info('Method: {}'.format(opts.dpx_method))
    # logger.info('FastQ: {}'.format(opts.input_fastq))
    logger.info('Barcode: {}'.format(opts.barcode_file))
    # logger.info('Target: {}'.format(opts.target_len))

    # logger.info('Misreadings_Barcode: {}'.format(opts.misreads_cutoff_cons))
    # logger.info('Misreadings_Constant: {}'.format(opts.misreads_cutoff_cons))
    # logger.info('Stats: {}'.format(opts.save_frequencies))

                                                                        ,opts.phredcutoff))
    try:
        os.makedirs(opts.data_path+'/Sequences/')
    except OSError:
        e = sys.exc_info()
        #print 'Warning, Folder',out_folder,' already exist'
        print 'Warning', e[1]

    for barcode in barcodes_list:
        trim(barcode, **opts.__dict__)

    return


if __name__ == '__main__':
    # Read argtments
    opts = get_options()

    # init logging
    time_stamp = time.ctime()
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s - \
                                %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename='run_{4}_{1}_{2}_{0}_{3}.log'.format(*time_stamp.split()),
                        filemode='w')
    logger = logging.getLogger(__name__)

    logger.info('JOB START {4} {1} {2} {0} {3}'.format(*time_stamp.split()))
    # DEMULTIPLEX
    workflow(opts)
    # DONE
    time_stamp = time.ctime()
    logger.info('JOB ENDS {4} {1} {2} {0} {3}'.format(*time_stamp.split()))
