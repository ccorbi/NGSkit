import os
import sys
import socket
try:
    hostname = socket.gethostname()
    if 'ccbr' in hostname or 'beagle' in hostname:
        sys.path.append('/home/kimlab2/ccorbi/python-Levenshtein-0.12.0/build/lib.linux-x86_64-2.7/')
except:
    pass

import Levenshtein as leven
#import argparse



def save_fastq(read1_id, read1_seq, read1_strand, read1_qual, l, suffix, out_dir):


    file_name = "split_{}_{}_F.fastq".format(suffix,l)
    out_file = os.path.join(out_dir, file_name).replace('\\', '/')
    f = open(out_file, 'a')
    f.write(read1_id)
    f.write(read1_seq)
    f.write(read1_strand)
    f.write(read1_qual)
    f.close()

    return


def is_longer_than(seq, l):

    extra = seq[l:]

    if extra[:10] == 'NNNNNNNNN':
        return False
    else:
        return True


def single_end(inputfile, out_dir='split_fastq', l=150,  **Kargs):
    '''.

        Quick method to split Fastq files base on reading length
    '''
    # Open Forward FastaQ file
    with open(inputfile, 'r') as read1:

            for read1_id in read1:
                # Read 4 by 4
                # ID lane info, seq info etc
                id1 = read1_id.partition(' ')
                # Read seq and Quality info
                sequence = ['', '']
                read1_seq = read1.next()
                read1_strand = read1.next()
                read1_qual = read1.next()
                # split base on
                # read1_id, read1_seq, read1_strand, read1_qual, l, suffix, out_dir
                if is_longer_than(read1_seq,l):
                    save_fastq(read1_id, read1_seq, read1_strand, read1_qual, l, 'long', out_dir)
                else:
                    save_fastq(read1_id, read1_seq, read1_strand, read1_qual, l, 'short', out_dir)


    return



if __name__ == '__main__':


    #opts = get_options()
    single_end(sys.argv[1],l=sys.argv[2])
