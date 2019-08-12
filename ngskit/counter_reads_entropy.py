import sys
import glob
import pandas as pd
import numpy as np
import argparse

from ngskit import analysis
from ngskit.utils import fasta_tools
from ngskit.utils import fastq_tools

from ngskit.utils import dna






def get_options():
    '''.

    '''
    parser = argparse.ArgumentParser(description="""
    Simple Reads counts


    Usage :
    %prog -f [Fasta file]  -f fastq

    """)


    parser.add_argument('-i', '--input_file', action="store", dest="input_f",
                        required=True, help='input fasta file from  demultiplex',
                        nargs='+')


    parser.add_argument('-o', '--output_file', action="store", dest="output_file",
                        required=True, help='Output name')

    parser.add_argument('-f', '--format', action="store", dest="iformat",
                        default="fasta",  help='Input file format, fasta, or fastq, (default: fastq)', type=str,
                        choices= ['fastq', 'fasta'])


    parser.add_argument('-e', '--encoding', action="store", dest="encoding",
                        default="NNS",  help='Input file format, fasta, or fastq, (default: fastq)', type=str,
                        choices= ['NNN', 'NNS','NNK'])


    
    parser.add_argument('-t', '--target-len', action="store", dest="target_len",
                        type=int, default=0, help='filter out sequences without the right len ')


    options = parser.parse_args()

    return options


def readfasta(input_fasta):

    dict_fasta = fasta_tools.read_fasta(input_fasta)

    return list(dict_fasta.values())



if __name__ == '__main__':

    # setup
    opts = get_options()

    

    fastqs = opts.input_f
    
    if opts.iformat== 'fasta':

        input_reader = readfasta
        cols = [ 'Seq']

    elif opts.iformat== 'fastq':

        input_reader = fastq_tools.read_fastq
        cols = ['Seq', 'avg_phred','phred']

    else:
        raise

    
    tdata = list()
    for fq in fastqs:
        a = input_reader(fq)
        tdata.extend(a)
        #print(a)

    df = pd.DataFrame(tdata, columns=cols)
    df.head()

    df['Ns'] = df['Seq'].str.contains('N')
    df['STOP'] = df['Seq'].apply(dna.code_4_any)
    # filter by previous criteria and translate to AA
    if opts.iformat== 'fastq':
        dfc = df[(df['STOP']==False)&(df['Ns']==False)&(df['avg_phred']>=30)]

        df_groupsdna =dfc.groupby('Seq',as_index=False).agg(len)
        del df_groupsdna['phred']
        del df_groupsdna['Ns']
        del df_groupsdna['STOP']
        df_groupsdna.rename(columns={'avg_phred':'Reads'}, inplace=True)

    
    if opts.iformat == 'fasta':
        dfc = df[(df['STOP']==False)&(df['Ns']==False)]

        df_groupsdna =dfc.groupby('Seq',as_index=False).agg(len)
        del df_groupsdna['STOP']
        df_groupsdna.rename(columns={'Ns':'Reads'}, inplace=True)



    df_groupsdna['aa'] = df_groupsdna['Seq'].apply(dna.translate2aa)

    df_groups_aa =df_groupsdna.groupby('aa')
    # Calc Shannon Entropy
    parsed_aa = list()
    for aa, df_grp in df_groups_aa:
        # if the group only have one member Check
        # if is variants = 1, otherwise, 0 Entropy
        

        if df_grp.shape[0] == 1:
            N = dna.possible_encondings(aa)
            if N == 1:
                # Seq, Reads, Entropy
                parsed_aa.append([aa, df_grp['Reads'].sum(),df_grp.shape[0], 1])
            else:
                parsed_aa.append([aa, df_grp['Reads'].sum(), df_grp.shape[0],0])
        else:
            N = dna.possible_encondings(aa)
            df_grp['prob'] = df_grp['Reads'] / df_grp['Reads'].sum(skipna=True)
            df_grp['eprob'] = df_grp['prob'] * (np.log2( df_grp['prob']))

            shannon_entropy = -1 * ( (1.0/(np.log2(N))) *  (df_grp['eprob'].sum(skipna=True)))

            # Seq, Reads, Entropy
            parsed_aa.append([aa, df_grp['Reads'].sum(),df_grp.shape[0], shannon_entropy])
    results =  pd.DataFrame(parsed_aa, columns=['Seq', 'Reads','Var' ,'E'])


    print('Total Sequences {}'.format(results.shape[0]))

    if opts.target_len:

        
        results['L'] = results['Seq'].str.len()
        results = results[results['L']==opts.target_len]

        print('Seqs After len filtering {}'.format(results.shape[0]))

    ent = results[results['E']>0]
    print('Seqs wiht E>0 {}'.format(ent.shape[0]))

    
    results.to_csv('{}_{}_preprocess.csv'.format(opts.output_file), index=False)
