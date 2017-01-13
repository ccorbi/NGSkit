from __future__ import print_function
import argparse
import operator

from ngskit.utils.dna import translate2aa

__author__ = "C. corbi-Verge"
__email__ = "carles.corbi@gmail.com"
__about__ = ''' This is a small script to count and rank reads from a
             fasta file. This script is agnostic, do not compare or map
             the reads with any source or library'''


def count_seqs_na(filename, seqs=dict()):
    '''.

    From a fasta file with na or aa sequences,
    count sequence and return a dict seq:counts

    '''
    with open(filename, 'r') as input_file:
        for line in input_file:
            if line.startswith('>'):
                pass
            else:
                seqs[line.strip()] = seqs.get(line.strip(), 0) + 1

    return seqs


def count_seqs_aa(filename, seqs=dict()):
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
                aa = translate2aa(line.strip())

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

    parser.add_argument('-n', help='Count oligos', dest='dna', default=False,
                        action='store_true')

    parser.add_argument('-o', '--output', action="store", dest="output",
                        required=True, help='Output Name')

    parser.add_argument('-f', '--input_fasta', action="store", dest="input_fasta",
                        required=True, help='input fasta file from  demultiplex',
                        nargs='+')

    options = parser.parse_args()

    return options


if __name__ == '__main__':

    opts = get_options()

    if opts.protein:
        counts = dict()
        for fname in opts.input_fasta:
            counts = count_seqs_aa(fname, counts)

        out_handelr = open('aa_' + opts.output, 'w')
        sorted_seqs = sorted(counts.items(), key=operator.itemgetter(1),
                             reverse=True)

        for seq in sorted_seqs:
            print('{}\t{}'.format(seq[0], seq[1]), file=out_handelr)

        out_handelr.close()

    if opts.dna:
        counts = dict()
        for fname in opts.input_fasta:
            counts = count_seqs_na(fname, counts)

        out_handelr = open('na_' + opts.output, 'w')

        sorted_seqs = sorted(counts.items(), key=operator.itemgetter(1),
                             reverse=True)

        for seq in sorted_seqs:
            print('{}\t{}'.format(seq[0], seq[1]), file=out_handelr)

        out_handelr.close()
