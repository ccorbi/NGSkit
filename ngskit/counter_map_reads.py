# -*- coding: utf-8 -*-

from __future__ import print_function
import pandas as pd
import argparse

from ngskit.utils.dna import translate2aa


__author__ = 'C. corbi-Verge'
__email__ = 'carles.corbi@gmail.com'
__about__ = \
    ''' This is a small script to count and rank reads from a fasta file. This script is agnostic, do not compare or map
             the reads with any source or library
bowtie-build -f WHEP.90k.rev-compl.108na.lib.fas WHEP.90k.rev-compl.108na.lib
for i in *CT*fasta;do echo ${i};bowtie -v 0 ../../../SEQ_PDZ_20170228/lib/cterm_7aa -f ${i} > ./bowtie/${i%.fasta}.out;done
counter_map_reads.py  --bowtie -o ${i%.out}.counter -i ${i} -s 50000;done

    '''


def load_blastout(filename):
    df = pd.read_csv(filename, delim_whitespace=True, names=[
        'queryId',
        'subjectId',
        'percIdentity',
        'alnLength',
        'mismatchCount',
        'gapOpenCount',
        'queryStart',
        'queryEnd',
        'subjectStart',
        'subjectEnd',
        'eVal',
        'bitScore',
        ])
    return df


def load_bowtieout(filename):
    df = pd.read_csv(filename, names=[
        'queryId',
        'strand',
        'subjectId',
        'off',
        'seq',
        'otheInstances',
        'quality',
        'mismatchDescriptor',
        ], sep='\t')

    # print(df.head(5))

    return df


def filter_blast(
    df,
    iden=100.0,
    l=54,
    m=0,
    ):
    df = df[(df['percIdentity'] >= 100.0) & (df['alnLength'] == l)
            & (df['mismatchCount'] <= m)]
    return df


def filter_bowtie(df, m=0):
    if m == 0:
        return df[df['mismatchDescriptor'].isnull()]


def return_cov(df, lib_size=6000):

    return df['subjectId'].unique().shape[0] / float(lib_size)


def from_bowtie_generate_ranking_file(df, outputname, **kargs):

    # Count Reads
    count = df.groupby(['subjectId'], as_index=False)['off'].agg(len)
    count.rename(columns={'off': 'Counts'}, inplace=True)
    # Count variants
    var = df.groupby(['subjectId' ], as_index=False)['off'].agg(len)
    var.rename(columns={'off': 'Var'}, inplace=True)
    # get the original seq and discard variants, I do not like it a lot
    pureseq = df[(df['off']==0) & (df['mismatchDescriptor'].isnull() )]
    pureseq = pureseq.drop_duplicates(subset=['subjectId'])
    pureseq = pureseq[['subjectId','seq']]
    # Merge
    final_data = pd.merge(pureseq, var, on= 'subjectId')
    final_data = pd.merge(final_data, count, on= 'subjectId')
    # sort and translate to AA
    final_data.sort_values(by='Counts', ascending=False, inplace=True)

    final_data['aa'] = final_data['seq'].apply(translate2aa)

    # print(q.head(5))

    final_data.to_csv(outputname, index=False, columns=['subjectId', 'aa',
             'Counts', 'Var'])
    return final_data


def get_options():

    parser = \
        argparse.ArgumentParser(description="""
    Simple Reads counts


    Usage Demultiplexation:
    %prog -f [Fasta file]  -p count AA space  -n count in NA space -o outputname

    """)

    parser.add_argument('--blast',
                        help='Count sequences from a blast output',
                        dest='blast', default=False, action='store_true'
                        )

    parser.add_argument('--bowtie',
                        help='Count Sequences from bowtie output',
                        dest='bowtie', default=False,
                        action='store_true')

    parser.add_argument(
        '-o',
        '--output',
        action='store',
        dest='output',
        required=True,
        help='Output Ranking file Name',
        )

    parser.add_argument(
        '-i',
        '--input_file',
        action='store',
        dest='input_file',
        required=True,
        help='file with mapped sequences bowtie/blast',
        )

    parser.add_argument(
        '-s',
        '--library_size',
        action='store',
        dest='library_size',
        help='Size of the library',
        type=int,
        )

    options = parser.parse_args()

    return options


if __name__ == '__main__':

    opts = get_options()

    if opts.blast:
        pass

    if opts.bowtie:

        data = load_bowtieout(opts.input_file)
        raw_seqs = data.shape[0]
        # this should be optional
        # usually I filter web i map to the library
        # data = filter_bowtie(data)
        filter_seqs = data.shape[0]
        unique_seqs = data['subjectId'].unique().shape[0]
        coverage = ''
        if opts.library_size:
            coverage = return_cov(data, lib_size=opts.library_size)

        if opts.output:
            g = from_bowtie_generate_ranking_file(data, opts.output)
            s = g['Counts'].describe()

            print('INPUT\tRaw Mapped READS\tFiltered Mapped READS\tUnique Reads\tlibrary Coverage\tmean\tstd\tmin\t25%\t50%\t75%\tmax'
                  )
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                opts.input_file,
                raw_seqs,
                filter_seqs,
                unique_seqs,
                coverage,
                s[u'mean'],
                s[u'std'],
                s[u'min'],
                s[u'25%'],
                s[u'50%'],
                s[u'75%'],
                s[u'max'],
                ))
        else:

            # print summary

            print('Raw Mapped READS\tFiltered Mapped READS\tUnique Reads\tlibrary Coverage\t'
                  )
            print('{}\t{}\t{}\t{}'.format(raw_seqs, filter_seqs,
                  unique_seqs, coverage))
