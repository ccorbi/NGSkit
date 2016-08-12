import os, sys
import socket

hostname = socket.gethostname()
if 'ccbr' in hostname or 'beagle'  in hostname:
    sys.path.append('/home/kimlab2/ccorbi/python-Levenshtein-0.12.0/build/lib.linux-x86_64-2.7/')

#from string import *
import Levenshtein as leven

import logging
import time
import pandas as pd
import optparse



#####
# Supplement



def match(seq, target, cutoff,read1_id,label=' '):

    '''Method used to compare sequence, such barcodes and constant regions.

    :param str seq: Nucleotide Sequence trimed from the fastaq file .
    :param str target: Template to compare, barcode or constant region.
    :param str cutoff: Maximum number of diff accepted between the two sequnces.
    :param str read1_id: Id of the sequences, to debug purposes. This info and the number of differences is recorder in the log file

    :return: False if there are diff over the cutoff
    :rtype: boolean'''


    cutoff = int(cutoff)
    distance = leven.distance(seq,target)

    if distance > cutoff:

        return False

    # Save to control purpuses
    if distance > cutoff and distance < 14:
        logger.info('Mismatch - {} - {} - {}'.format(label,read1_id.strip(),distance))

    return True

def read_barcodes(barcode_file):

    '''Method to Read barcode file, all the seq must be in 5' 3' sense

    :param str barcode_file: The path or name of the file that contains barcode info

    :return: a dict of barcode, name is the sample_id for a list of sequences.
    :rtype: dict '''

    assert os.path.isfile(barcode_file)

    barcode = {}

    with open(barcode_file,'r') as input_file:
        for line in input_file:
            line = line.strip()

            if not line.startswith('#'):
                data = line.split()
                if len(data) == 5:
                    name = data[0].strip()
                    barcode[name] = [data[1].strip(), data[2].strip(),data[3].strip(), data[4].strip()] #Forward-Barcode, Barcode-2-Reversed, Konstant_region1, Konstant_region2-Reversed
                else:
                    print 'Barcode should contain #sample ID, #Forward-Barcode, Barcode-2-Reversed, Konstant_region1, Konstant_region2-Reversed'
                    print  data

    return barcode

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
    #Open Forward
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

def counting_reads(blast_output,aln_len=0,filter_by_lib='',quality_check=False):

    ''''Read the blast output, and count reads that comes from the blast that are in the library, perfect match [need improvment]

    blastall -p blastn -d design_lib.fasta -i smaples_id_demultiplex.fasta -o blastoutput -v1 -b1 -m8

    /home/kimlab2/ccorbi/ncbi-blast/bin/blastn -db /home/kimlab2/ccorbi/optim_lib/Kim/fasta_libs/optm_library.fasta -query ${i}  -evalue 1 -out ${i}.out  -outfmt '6 qacc sacc pident length  mismatch gapopen qstart qend sstart send evalue bitscore'

    :param str blastoutput: Path to blast output
    :param str aln_len: length of the oligos in the library, if is not provided, the length is not check

    :param str filter_by_lib: This is  a special function, do it ad hoc, for one time I have to screen without the constant regions.

    :return: print Counts
    '''

    #blastoutput = param.get('blast_output','')

    qlty_landscape = {}

    blast_info = pd.read_csv(blast_output,delim_whitespace=True,names=['queryId', 'subjectId', 'percIdentity', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore'])

    qlty_landscape['percIdentity'] = blast_info['percIdentity']
    qlty_landscape['mismatchCount'] = blast_info['mismatchCount']

    if aln_len==0:
        good_data =   blast_info[blast_info['percIdentity']==100.0]
    else:
        good_data =   blast_info[(blast_info['percIdentity']==100.0)&(blast_info['alnLength']==aln_len)]

    if filter_by_lib!='':
        #read fasta, and filter by keys
        lib = read_fasta2dict(filter_by_lib)
        good_data = good_data[good_data['subjectId'].isin(lib.keys())]


    if quality_check:
        print qlty_landscape['percIdentity'].value_counts()
        print qlty_landscape['mismatchCount'].value_counts()

    # print 'File\tTotal hits\tUnique hits'
    print blast_output,good_data.shape[0],len(good_data['subjectId'].unique())

    return


class NGS_jobs(object):

    """docstring for Demultiplex. Main class to deal with peptides screening data [need improvement]
    """


    def __init__(self, **param):

        ''' Init class
        :param str inputfile : fastq data with the seq
        :param str out_path: Name of the output path,
                              demultiplex_output (default)
        :param str design_len: length of the design sequences
        :param str phredcutoff: Reading Quality cutoff (default 30)
        :param str cutoff_cons: Max number of misreading allowed in
                               the constant constant_region (default 1)
        :param str cutoff_barcode: Max number of misreading allowed in
                                  the constant constant_region  (default 1) '''

        self.inputfile = param.get('inputfile','')
        self.peplength = param.get('design_len',0)
        self.out_path = param.get('out_path','demultiplex_output')

        # self.demultiplexed_file = param.get('demultiplexed_file','')

        self.reading_mismatch_cutoff_cons = param.get('cutoff_cons', 1)
        self.reading_mismatch_cutoff_barcode = param.get('cutoff_barcode', 1)
#
#        self.barcodes = param.get('barcodes', {})
        # self.sample_id = param.get('sample_id', {})

        self.reading_quality_cutoff = param.get('phredcutoff', 30)
        #Counting Sequence paramater
        self.quality_check = param.get('quality_check', False)



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


    def demultpx_(self,inputfile,barcode_file,out_path='demultiplex',dpx_type='QUICK',dump=False,**Kargs):
        '''Main method to demultiplex Single end sequences from deepSeq. No quality score applyed [need improvements]

        :param str inputfile: Fastaq file
        :param str barcode_file: Path to the file with the barcodes and constant regions
        :param str out_path: path to save demultiplexed fastaq files

        :param str dpx_type: Type of demultiplexation by default quick
                QUICK: Only the first barcode and constant region are checked
                STANDARD: Both barcodes and constant regions are check
                BARCODES: Only the barcodes are used
                FS: frame shift search, Flexible search of the second constant region and barcode'''

        #Check if the type is right
        if  not dpx_type in ['QUICK', 'STANDARD', 'BARCODES','FS']:
            print dpx_type,'do not match with any type of demultiplex'
            print 'Accepted types : QUICK, STANDARD, BARCODES, FS (frame_shift)'
            print 'QUICK is used by default'
            dpx_type='QUICK'


        #load barcodes file

        barcodes = read_barcodes(barcode_file)

        #make output folder

        makeoutputdirs(barcodes,out_path)

        #check barcodes integrity, peplength, fastq
        assert len(barcodes.items())>0
        assert int(self.peplength) >0
        #print in the log

        #Check if something must be dump

        if dump:
            o = './dump/{}/'.format(out_path)
            try:
                os.makedirs(o)
            except OSError:
                e = sys.exc_info()
            #print 'Warning, Folder',out_folder,' already exist'
                print 'Warning', e[1]

            self.dump_d = dict()
            self.dump_cons1 =dict()
            self.dump_cons2 = dict()
            for sample_id in barcodes.keys():
                self.dump_cons1[sample_id]  = open(o+sample_id+'_K1.dump.out','w')
                self.dump_cons2[sample_id]  = open(o+sample_id+'_K2.dump.out','w')
                self.dump_d[sample_id]  = open(o+sample_id+'_D.dump.out','w')

        #if we use only barcodes to demultiplex
        #is like a statd but we do not check the constant region identity.
        if dpx_type=='BARCODES':
            self.reading_mismatch_cutoff_cons = 100000
            dpx_type='STANDARD'

        # #Check barcode


        #Init File Stats kind of log
        self.stats_init(inputfile,barcodes)




        with open(inputfile, 'r') as read1:
            #Open Reverse
            # with open(read2_path, 'r') as read2:
                for read1_id in read1:
                    id1 = read1_id.partition(' ')

                    sequence = ['', '']
                    read1_seq = read1.next()
                    read1.next()
                    read1_qual = read1.next()


                    for sample_id in barcodes:
                        barcode1 = barcodes[sample_id][0]
                        barcode2 = barcodes[sample_id][1]


                        F_cons = barcodes[sample_id][2]
                        S_cons = barcodes[sample_id][3]

                        bar1_len = len(barcode1)
                        bar2_len = len(barcode2)

                        cons1_len = len(F_cons)
                        cons2_len = len(S_cons)



                        #extract barcodes from seq and compare
                        # B1 + K1 + pep + K2 + B2
                        #save when seq match in a different length

                        #Extract Bacode 1, always index = 0
                        b1 = read1_seq[0:bar1_len]
                        #Extract constant region 1, next to barcode 1
                        cons1 = read1_seq[bar1_len:cons1_len+bar1_len]

                        #cons2 = read1_seq[len(barcode1)+len(F_cons)+int(peplength):len(barcode1)+len(F_cons)+int(peplength)+len(S_cons)]
                        if match(b1, barcode1, self.reading_mismatch_cutoff_barcode ,read1_id,label='b1'):
                            #file stats
                            self.stats_write(inputfile,sample_id,barcode1)


                            if match(cons1, F_cons, self.reading_mismatch_cutoff_cons,read1_id,label='c1'):
                                #Floating window ... some seq can be shorter, value info
                                #first find the Second constannat region
                                self.stats_write(inputfile,sample_id,F_cons)
                                #
                                if dump:
                                    print>>self.dump_cons1[sample_id],cons1.strip()


                                if dpx_type=='QUICK':
                                        self.stats_write(inputfile,sample_id,'Hit')
                                        self.save_seq(read1_id,read1_seq,read1_qual,sample_id,self.peplength,out_path)
                                        if dump:
                                            design = read1_seq[bar1_len+cons1_len:bar1_len+cons1_len+int(self.peplength)]
                                            cons2 = read1_seq[bar1_len+cons1_len:bar1_len+cons1_len+int(self.peplength)+cons2_len]
                                            print>>self.dump_d[sample_id],design
                                            print>>self.dump_cons2[sample_id],cons2





                                #frame shift
                                elif dpx_type=='FS':
                                    # cons2 = read1_seq[bar_len1+cons_len1+int(self.peplength):bar_len1+cons_len1+int(self.peplength)+cons_len2]

                                    fshift = self.look_4_frameshift(read1_seq,read1_id,bar1_len,cons1_len,cons2_len,S_cons,bar2_len,barcode2,sample_id,inputfile,dump)

                                    if fshift != None:
                                        #Save Stats
                                        self.stats_write(inputfile,sample_id,'Hit')
                                        self.stats_write(inputfile,sample_id,fshift)

                                        #Save Sequence
                                        fshift = int(self.peplength) + fshift
                                        self.save_seq(read1_id,read1_seq,read1_qual,sample_id,fshift,out_path)



                                elif dpx_type=='STANDARD':

                                    cons2 = read1_seq[bar1_len+cons1_len+int(self.peplength):bar1_len+cons1_len+int(self.peplength)+cons2_len]
                                    if match(cons2, S_cons, self.reading_mismatch_cutoff_cons,read1_id,label='c2'):


                                        #Save Stadistics
                                        self.stats_write(inputfile,sample_id,S_cons)
                                        #Once Constant region it is located, extract barcode 2, using shift variable s

                                        if dump:

                                            design = read1_seq[bar1_len+cons1_len:bar1_len+cons1_len+int(self.peplength)]
                                            #cons2 = read1_seq[bar1_len+cons1_len:bar1_len+cons1_len+int(self.peplength)+cons2_len]
                                            print>>self.dump_d[sample_id],design
                                            print>>self.dump_cons2[sample_id],cons2

                                        b2 = read1_seq[bar1_len+cons1_len+int(self.peplength)+cons2_len:bar1_len+cons1_len+int(self.peplength)+cons2_len+bar2_len]
                                        if match(b2, barcode2, self.reading_mismatch_cutoff_barcode,read1_id,label='b2'):

                                            self.stats_write(inputfile,sample_id,barcode2)
                                            self.stats_write(inputfile,sample_id,'Hit')

                                            self.save_seq(read1_id,read1_seq,read1_qual,sample_id,self.peplength,out_path)





        #write output.

        self.gbl_stats.to_csv(out_path+'-'+barcodes.keys()[0]+'_gbl_stats.csv')

        if dump:
            for segment in [self.dump_cons1,self.dump_cons2,self.dump_d]:
                for dump_file in  segment.itervalues():
                    dump_file.close()



    def look_4_frameshift(self,read1_seq,read1_id,bar1_len,cons1_len,cons2_len,S_cons,bar2_len,barcode2,sample_id,inputfile,dump):


        retro_start = -10

        for i in range(26):
            shift = retro_start + i

            #extract constant region 2, using shift variable to scan 10 positions
            cons2 = read1_seq[bar1_len+cons1_len+int(self.peplength)+shift:bar1_len+cons1_len+int(self.peplength)+shift+cons2_len]

            if match(cons2, S_cons, self.reading_mismatch_cutoff_cons,read1_id,label='c2'):
                frame_shift = int(self.peplength) + shift


                #Save Stadistics
                self.stats_write(inputfile,sample_id,S_cons)
                #Once Constant region it is located, extract barcode 2, using shift variable s

                if dump:

                    design = read1_seq[bar1_len+cons1_len:bar1_len+cons1_len+int(self.peplength)+shift]

                    print>>self.dump_d[sample_id],design
                    print>>self.dump_cons2[sample_id],cons2

                b2 = read1_seq[bar1_len+cons1_len+cons2_len+frame_shift:bar1_len+cons1_len+cons2_len+frame_shift+bar2_len]
                if match(b2, barcode2, self.reading_mismatch_cutoff_barcode,read1_id,label='b2'):
                    #Save Stadistics
                    self.stats_write(inputfile,sample_id,barcode2)


                    return shift

        return None

    def save_seq(self,read1_id,read1_seq,read1_qual,sample_id,frame_shift,out_path):


        file_name = sample_id + "_"+str(frame_shift)+"_F.fastq"
        out_file = os.path.join(out_path, sample_id, file_name).replace('\\', '/')
        f = open(out_file, 'a')
        f.write(read1_id)
        f.write(read1_seq)
        f.write("+\n")
        f.write(read1_qual)
        f.close()

        return




    def stats_init(self,inputfile,barcode):

        '''
        Format Dataframe to save Reads

        '''

        total_ids = []
        for sample_id in barcode:

            total_ids.append(barcode[sample_id][0])
            total_ids.append(barcode[sample_id][1])
            total_ids.append(barcode[sample_id][2])
            total_ids.append(barcode[sample_id][3])

        # Barcodes Constant Regions + FrameShifts
        columns_id = list(set(total_ids))+['-10','-9','-8','-7','-6','-5','-4',
        '-3','-2','-1','0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']

        #init stats
        columns_id.append('Hit')
        rows_id = barcode.keys()
        rows_id.append(inputfile)

        gbl_stats = pd.DataFrame(index=rows_id,columns=columns_id)
        gbl_stats.fillna(0,inplace=True)

        self.gbl_stats = gbl_stats

        return

    def stats_write(self,inputfile,sample_id,segment):
        '''Tag counts '''
        self.gbl_stats.at[inputfile,str(segment)] = self.gbl_stats.at[inputfile,str(segment)] + 1
        self.gbl_stats.at[sample_id,str(segment)] = self.gbl_stats.at[sample_id,str(segment)] +1

        return

    def dump(self,option):
        '''This method dump data to a file '''


        return




def get_options():

    parser = optparse.OptionParser()


    parser.set_usage("""
    Deep Sequencing extraction sequences tools


    Usage Demultiplexation:
    %prog -a demultpx -b [BarCode_file.inp]  -i [deep_seq_file.fastq] -o [folder_name] -l 54 -t QUICK -options

    Usage Extraction Sequences:
    %prog -a extract -b [BarCode_file.inp]  -d [data_folder]  -l 54 -options
    """)

    parser.add_option('-a', '--action',
                      type='choice',
                      action='store',
                      dest='action',
                      choices=['demultpx', 'extract', 'testing','counting'],
                      default=False,
                      help='What you want to do: demultpx, extract_seq, testing')

    parser.add_option('-l','--design_len', action="store", dest="design_len",
        default=False,help='Lengh of the designed oligo')

    parser.add_option('-i','--input_file', action="store", dest="inputfile",
        default=False,help='inputfile FASTAQ (demultiplex) or blastoutput (counting)')

    parser.add_option('-b','--barcode_file', action="store", dest="barcode_file",
        default=False,help='File that contains barcodes and cosntant regions')

    parser.add_option('-o','--out_path', action="store", dest="out_path",default='demultiplex',
        help='Output folder, called demultiplex by default')


    parser.add_option('-d','--data_path', action="store", dest="data_path",default='demultiplex',
        help='Data folder, called demultiplex by default where the FASTAQ demultiplxed files are saved. Only need it for the Extraction action')

    # optional
    parser.add_option('-t','--demultiplex_type', action="store", dest="dpx_type",default='STANDARD',
        type='choice',choices=['QUICK', 'STANDARD', 'BARCODES','FS'],
        help="""Type of demultiplexation by default quick; \
        QUICK, Only the first barcode and constant region are checked \
        STANDARD, Both barcodes and constant regions are check\
        BARCODES, Only the barcodes are used \
        FS: frame shift search, Flexible search of the second constant region and barcode""")

    parser.add_option('--cutoff_cons', action="store", dest="cutoff_cons", default=1,type='int',
        help='Max number of misreading allowed in the constant constant_region (default 1)')

    parser.add_option('--cutoff_barcode', action="store", dest="cutoff_barcode", default=1,type='int',
        help='Max number of misreading allowed in the constant constant_region  (default 1)')

    parser.add_option('--dump', help='Dump constant regions', dest = 'dump', default = False,
                      action = 'store_true')

    parser.add_option('--cutoff_quality', help='Cutoff Quality of the fastq, default 30 ',
                      default = 30, dest="phredcutoff",type='int', action = 'store')

    parser.add_option('-c','--chopp', action="store", dest="chopp",default='peptide',
        type='choice',choices=['peptide', 'constant'], help='When extract seq, only peptide region or plus constant regions')

    (options, args) = parser.parse_args()




    # Check dependencies
    if options.action is False:   # if action is not given
        parser.print_help()
        parser.error('Action not given')


    if options.action == 'demultpx':
        if options.design_len and options.inputfile and options.barcode_file:
            pass
        else:
            parser.print_help()
            print ''
            parser.error('For demultiplexation the following fields are mandatory; Sequencing file (-i), barcode file (-b) and Design length (-l)')

    if options.action == 'extract':
        if options.design_len and options.data_path and options.barcode_file:
            pass
        else:
            parser.print_help()
            print ''
            parser.error('For Extraction the following fields are mandatory; Data Folder (-d), barcode file (-b) and Design length (-l)')


    return options

def main(opts):



    if opts.action == 'demultpx':


#        job =  NGS_jobs(design_len=opts.design_len)

        job =  NGS_jobs(**opts.__dict__)
        # Here Pass information from options
        logger.info('Start {} with method {} on {};barcode file {};Design length {}'.format(opts.action
                                                                                    ,opts.dpx_type
                                                                                    ,opts.inputfile
                                                                                    ,opts.barcode_file
                                                                                    ,opts.design_len))
        # Call the action
        job.demultpx_(**opts.__dict__)

    if opts.action == 'extract':
        logger.info('Start {}  on {};barcode file {};Design length {} Chopp {}  Q {}'.format(opts.action
                                                                            ,opts.data_path
                                                                            ,opts.barcode_file
                                                                            ,opts.design_len
                                                                            ,opts.chopp
                                                                            ,opts.phredcutoff))
        try:
            os.makedirs(opts.data_path+'/Sequences/')
        except OSError:
            e = sys.exc_info()
            #print 'Warning, Folder',out_folder,' already exist'
            print 'Warning', e[1]
        ### need to improve this
        job =  NGS_jobs(**opts.__dict__)
        job.trim_all_(**opts.__dict__)

    if opts.action == 'counting':

       logger.info('start {} on '.format(opts.action,opts.inputfile))
       counting_reads(blast_output=opts.inputfile)

    return



if __name__ == '__main__':


    opts = get_options()

    #init logging
    logger = logging.getLogger('{}'.format(opts.action))
    logger.setLevel(logging.DEBUG)
    # # create file handler which logs even debug messages
    # Collect them in the same directory is much better solution to control them
    start_time = time.time()
    fh = logging.FileHandler('run_'+str(start_time)+'.log')
    fh.setLevel(logging.DEBUG)
    # # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    main(opts)