from __future__ import print_function
import argparse
import operator

from ngskit.utils.dna import translate2aa
from ngskit.utils.fasta_tools import write_fasta_sequence

__author__ = "C. corbi-Verge"
__email__ = "carles.corbi@gmail.com"
__about__ = ''' This is a small script to count and rank reads from a
             fasta file. This script is agnostic, do not compare or map
             the reads with any source or library'''


def count_seqs_na(filename, seqs=dict(), min_lenght=0):
    '''.

    From a fasta file with na or aa sequences,
    count sequence and return a dict seq:counts

    '''
    with open(filename, 'r') as input_file:
        for line in input_file:
            if line.startswith('>'):
                pass
            else:
                # do not count seq under a length, usefull to filter framshifts, etc
                if min_lenght and min_lenght <= len(line.strip()):
                    pass
                else:
                    seqs[line.strip()] = seqs.get(line.strip(), 0) + 1

    return seqs


def count_seqs_aa(filename, seqs=dict(), min_lenght=0):
    '''.

     From a fasta file with na seq,
    transform to aminoacids and count sequences,
    finally  return a dict seq:counts

    '''
    with open(filename, 'r') as input_file:
        for line in input_file:
            if line.startswith('>'):
                pass
            else:
                try:
                    aa = translate2aa(line.strip())

                except KeyError:
                    # so far ignore it found a N
                    # on the future add fuzzy logic to add read to a sequence
                    continue
                # do not count seq under a length, usefull to filter framshifts, etc
                if min_lenght and min_lenght <= len(aa):
                    pass
                else:
                    seqs[aa] = seqs.get(aa, 0) + 1
                # print >> fasta_2_logo,translate2aa(line.strip())
    # sorted_seqs = rank_reads(seqs)
    return seqs


def get_options():
    '''.

    '''
    parser = argparse.ArgumentParser(description="""
    Simple Reads counts


    Usage Demultiplexation:
    %prog -f [Fasta file]  -p count AA space  -n count in NA space -o outputname

    """)

    parser.add_argument('-p', help='Translate Oligos to Aminoacids and count',
                        dest='protein', default=False, action='store_true')

    parser.add_argument('-o', '--output', action="store", dest="output",
                        required=True, help='Output Name')

    parser.add_argument('-i', '--input_fasta', action="store", dest="input_fasta",
                        required=True, help='input fasta file from  demultiplex',
                        nargs='+')

    parser.add_argument('-f', '--output_format', action="store", dest="output_format",
                        default="csv",  help='Output as  "fasta",by default  "csv" file', type=str,
                        choices= ['csv', 'fasta'])

    parser.add_argument('-m', '--min-len', action="store", dest="min_lenght",
                        type=int, default=0, help='Ignore sequence under this lenghth')

    options = parser.parse_args()

    return options


def write_output(sequences, fileouthandel,  otype ):

    for seq in sequences:

        if otype == "csv":
            print('{}\t{}'.format(seq[0], seq[1]), file=fileouthandel)
        elif otype == "fasta":
            write_fasta_sequence([ seq[0]+'_'+str(seq[1]), seq[0]], fileouthandel)

    fileouthandel.close()


if __name__ == '__main__':

    # setup
    opts = get_options()

    if opts.protein:
        counter_fn = count_seqs_aa
        output_prefix = 'aa_'

    else:
        counter_fn = count_seqs_na
        output_prefix = 'na_'

    # do the job for each file
    counts = dict()
    for fname in opts.input_fasta:
        counts = counter_fn(fname, counts, min_lenght=opts.min_lenght)

    out_handelr = open(output_prefix + opts.output, 'w')
    sorted_seqs = sorted(counts.items(), key=operator.itemgetter(1),
                         reverse=True)

    write_output(sequences=sorted_seqs, fileouthandel=out_handelr, otype=opts.output_format)
