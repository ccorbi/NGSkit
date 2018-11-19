# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division
import sys
import glob
import pandas as pd
import numpy as np

from  ngskit.utils import dna
#form dna_util import *
from common import *

# Pipelines
def lentivirus_combine(data_path = '/home/ccorbi/Work/Beagle/optim_lib/Kim/demultpx_barcodeOnly/Sequences/',
               gaps='fill', fill_value=0, time_points=[0,3,7,14],
                cell_times = {'HEK293T':[0,1.8,3.9,9.4],'RWP1':[0,2.0,4.1,9.7]}):
    """Read data from the blast ouput files. Apply DESEQ normalization.
    Merge Replicas and average the reads, calc dropout score.
    One library at a time.

    :param str data_path: path with the where the files are located
    :param str gaps: Behaivor with the gaps, by default fill, fill the gaps with one read.
    :param func normfunction: Function to normalize the reads, by default log10
    :param list time_points: list contain the time points info [0,1,4,10...]
    :param dict cell_times: Dictionary where key is the cell lines name and
    the values a list with the Double time of the cell line for each time point. [0,2.4,4.5,10]
    This must match the time points.

    :param bool normed: Divison all the time points by the biggest value in the time serie
    :param bool itemization: return slope for each time point and the dropscore
    :param int scalar: scalar to multiply dropscore [increse the signal]

    :returns:
        :rtype dict: The information for each lib and cell lines is merge toghether in a
        dataframe. Each DataFrame is placed in a dictionary with key [lib_cellline]
    """
    alldata = dict()


    # Check if the number of times match the doubling times
    for times in cell_times.values():
        assert len(times) == len(time_points)


    for cellline in cell_times.keys():
        print('Reading {}'.format(cellline))

        #get reads for all replicas in time point
        # init vars
        # DataFrame
        poll_of_replicas = pd.DataFrame(columns=['referenceId'])
        # recollect columns with reads info, for normalize
        columns_reads = list()
        # for group replicas (average)
        columns_replicas = dict()

        for t in time_points:
            # init dict

            timereplicas = glob.glob(data_path+'*%s*T%i*.out' % (cellline,t))
            columns_replicas[t] = list()

            #Get reads for replica, normal and put it toghet
            # Read all replicas for a T point
            for replica in timereplicas:
                print('Parsing Time {} - replica {}'.format(t, replica))

                # Read data for replica
                # No simple norm
                raw_data = read_bowtie(replica, norm=False)

                # Get name replica and label it
                name_replica = replica.split('/')[-1]
                columns_reads.append('Reads_{}'.format(name_replica))
                columns_replicas[t].append('nReads_{}'.format(name_replica))
                raw_data.rename(columns={'Reads':'Reads_{}'.format(name_replica)}, inplace=True)
                #                         'nReads':'nReads_{}'.format(name_replica)}, inplace=True)
                # Join replicas

                poll_of_replicas = pd.merge(poll_of_replicas, raw_data, on=['referenceId'], how='outer')


        # Handle the Gaps in the Time points Not to the replicas.
        if gaps == 'fill':
            poll_of_replicas.fillna(fill_value, inplace=True)
        # if gaps is not set up
        # Drop row where there is some missing values
        else:
            poll_of_replicas.dropna(inplace=True)


        # apply DESEQ normalization
        # consider all samples like replicas
        print(poll_of_replicas.shape[0])
        poll_of_replicas = norm_deseq(poll_of_replicas, columns_reads)
        print(poll_of_replicas.shape[0])

        # join & Average all replicas
        # this should be done after the normalization
        for t in time_points:
            poll_of_replicas['nReads_T{}'.format(t)] = poll_of_replicas[columns_replicas[t]].mean(axis=1, skipna=True)



        # calc dropOut score log
        assert(poll_of_replicas.shape[0] > 0)

        poll_of_replicas['DPlog_{}'.format(cellline)] = poll_of_replicas.apply(dropOutScoreLog2,
                                                                            sampling_points=time_points,
                                                                            dbl_times = cell_times[cellline],
                                                                            axis=1)
        # Calc Lienal regression
        poll_of_replicas['LR_{}'.format(cellline)] = poll_of_replicas.apply(DP_linreg,
                                                                    sampling_points=time_points,
                                                                    dbl_times = cell_times[cellline],
                                                                    axis=1)


        alldata[cellline] = poll_of_replicas.copy()
        # print some values to see it is working properlly
        print('DONE {}'.format(cellline))
        print('{}     '.format(alldata[cellline].shape[0]))
        print(alldata[cellline].head(5))


    return alldata


def lentivirus_deseq(data_path = '/home/ccorbi/Work/Beagle/optim_lib/Kim/demultpx_barcodeOnly/Sequences/',
               gaps='fill', fill_value=0, time_points=[0,3,7,14],
                cell_times = {'HEK293T':[0,1.8,3.9,9.4],'RWP1':[0,2.0,4.1,9.7]}):
    """Read data from the blast ouput files. Apply DESEQ normalization.
    Merge Replicas and average the reads, calc dropout score.
    One library at a time.

    :param str data_path: path with the where the files are located
    :param str gaps: Behaivor with the gaps, by default fill, fill the gaps with one read.
    :param func normfunction: Function to normalize the reads, by default log10
    :param list time_points: list contain the time points info [0,1,4,10...]
    :param dict cell_times: Dictionary where key is the cell lines name and
    the values a list with the Double time of the cell line for each time point. [0,2.4,4.5,10]
    This must match the time points.

    :param bool normed: Divison all the time points by the biggest value in the time serie
    :param bool itemization: return slope for each time point and the dropscore
    :param int scalar: scalar to multiply dropscore [increse the signal]

    :returns:
        :rtype dict: The information for each lib and cell lines is merge toghether in a
        dataframe. Each DataFrame is placed in a dictionary with key [lib_cellline]
    """
    alldata = dict()


    # Check if the number of times match the doubling times
    for times in cell_times.values():
        assert len(times) == len(time_points)


    for cellline in cell_times.keys():
        print('Reading {}'.format(cellline))

        #get reads for all replicas in time point
        # init vars
        # DataFrame
        poll_of_replicas = pd.DataFrame(columns=['referenceId'])
        # recollect columns with reads info, for normalize
        columns_reads = list()
        # for group replicas (average)
        columns_replicas = dict()

        for t in time_points:
            # init dict

            timereplicas = glob.glob(data_path+'*%s*T%i*.out' % (cellline,t))
            columns_replicas[t] = list()

            #Get reads for replica, normal and put it toghet
            # Read all replicas for a T point
            for replica in timereplicas:
                print('Parsing Time {} - replica {}'.format(t, replica))

                # Read data for replica
                # No simple norm
                raw_data = read_bowtie(replica, norm=False)

                # Get name replica and label it
                name_replica = replica.split('/')[-1]
                columns_reads.append('Reads_{}'.format(name_replica))
                columns_replicas[t].append('nReads_{}'.format(name_replica))
                raw_data.rename(columns={'Reads':'Reads_{}'.format(name_replica)}, inplace=True)
                #                         'nReads':'nReads_{}'.format(name_replica)}, inplace=True)
                # Join replicas

                poll_of_replicas = pd.merge(poll_of_replicas, raw_data, on=['referenceId'], how='outer')


        # Handle the Gaps in the Time points Not to the replicas.
        if gaps == 'fill':
            poll_of_replicas.fillna(fill_value, inplace=True)
        # if gaps is not set up
        # Drop row where there is some missing values
        else:
            poll_of_replicas.dropna(inplace=True)


        # apply DESEQ normalization
        # consider all samples like replicas
        print(poll_of_replicas.shape[0])
        poll_of_replicas = norm_deseq(poll_of_replicas, columns_reads)
        print(poll_of_replicas.shape[0])

        # join & Average all replicas
        # this should be done after the normalization
        for t in time_points:
            poll_of_replicas['nReads_T{}'.format(t)] = poll_of_replicas[columns_replicas[t]].mean(axis=1, skipna=True)



        # calc dropOut score
        assert(poll_of_replicas.shape[0] > 0)

        poll_of_replicas['DS_{}'.format(cellline)] = poll_of_replicas.apply(dropOutScore,
                                                                            sampling_points=time_points,
                                                                            dbl_times = cell_times[cellline],
                                                                            axis=1)


        alldata[cellline] = poll_of_replicas.copy()
        # print some values to see it is working properlly
        print('DONE {}'.format(cellline))
        print('{}     '.format(alldata[cellline].shape[0]))
        print(alldata[cellline].head(5))


    return alldata

def lentivirus(data_path = '/home/ccorbi/Work/Beagle/optim_lib/Kim/demultpx_barcodeOnly/Sequences/',
               gaps='fill', normfunction=np.log10, fill_value=0,
               libs = ['Cterminal','Disorderome'],time_points=[0,3,7,14],
                cell_times = {'HEK293T':[0,1.8,3.9,9.4],'RWP1':[0,2.0,4.1,9.7]}):
    """ read data from the blast ouput files. Merge Replicas and Normalize the reads
    (by default: Log10 of the reads, and avg of the replicas). collect different
    time point and build a DataFrame with the info.

    :param str data_path: path with the where the files are located
    :param str gaps: Behaivor with the gaps, by default fill, fill the gaps with one read.
    :param func normfunction: Function to normalize the reads, by default log10
    :param list time_points: list contain the time points info [0,1,4,10...]
    :param dict cell_times: Dictionary where key is the cell lines name and
    the values a list with the Double time of the cell line for each time point. [0,2.4,4.5,10]
    This must match the time points.

    :param bool normed: Divison all the time points by the biggest value in the time serie
    :param bool itemization: return slope for each time point and the dropscore
    :param int scalar: scalar to multiply dropscore [increse the signal]

    :returns:
        :rtype dict: The information for each lib and cell lines is merge toghether in a
        dataframe. Each DataFrame is placed in a dictionary with key [lib_cellline]
    """
    alldata = dict()


    # Check if the number of times match the doubling times
    for times in cell_times.values():
        assert len(times) == len(time_points)

    for lib in libs:
        for cellline in cell_times.keys():
            print('Reading {} {}'.format(lib, cellline))

            #get reads for all replicas in time point
            # init vars
            poll_of_replicas = pd.DataFrame(columns=['referenceId'])
            columns_replicas = dict()

            for t in time_points:
                # init dict

                timereplicas = glob.glob(data_path+'*%s*T%i*.out' % (cellline,t))
                columns_replicas[t] = list()

                #Get reads for replica, normal and put it toghet
                # Read all replicas for a T point
                for replica in timereplicas:
                    print('Parsing Time {} - replica {}'.format(t, replica))

                    # Read data for replica
                    # In this library, one NA seq per peptide, so we can normilize here
                    raw_data = read_bowtie(replica, norm=norm_TC)
                    # Get name replica
                    name_replica = replica.split('/')[-1]
                    columns_replicas[t].append('nReads_{}'.format(name_replica))
                    raw_data.rename(columns={'Reads':'Reads_{}'.format(name_replica),
                                             'nReads':'nReads_{}'.format(name_replica)}, inplace=True)
                    # Join replicas
                    poll_of_replicas = pd.merge(poll_of_replicas, raw_data, on=['referenceId'], how='outer')


                # join & Average all replicas
                poll_of_replicas['nReads_T{}'.format(t)] = poll_of_replicas[columns_replicas[t]].mean(axis=1,skipna=True)

            # Handle the Gaps in the Time points Not to the replicas.
            if gaps == 'fill':
                poll_of_replicas.fillna(fill_value,inplace=True)
            # Drop row where there is some missing values
            else:
                poll_of_replicas.dropna(inplace=True)


            # calc dropOut score
            assert(poll_of_replicas.shape[0] > 0)

            poll_of_replicas['DS_{}'.format(cellline)] = poll_of_replicas.apply(dropOutScore,
                                                                                sampling_points=time_points,
                                                                                dbl_times = cell_times[cellline],
                                                                                axis=1)


            alldata[lib+'_'+cellline] = poll_of_replicas.copy()
            # print some values to see it is working properlly
            print('DONE {}'.format(lib+'_'+cellline))
            print('{}     '.format(alldata[lib+'_'+cellline].shape[0]))
            print(alldata[lib+'_'+cellline].head(5))


    return alldata


def b2h_random(fasta_file, norm=True):
    """Pipeline for preprocess of random libraries from b2h.

        Read demultiplexed fasta files, filter non NNS and Frameshift sequences.
    Returns a Dataframe with the sequences group by AA

    Parameters
    ----------
        fasta_file : str
            Path to the fasta file

    Returns
    -------
        Dataframe
            Return a dataframe with filtered sequneces
            ['referenceId', 'Seq', 'Reads', 'nReads']


    """
    # load data
    df = read_fasta(fasta_file, norm=False, translate2aa=False)
    # remove extra stop codons
    # df['na'] = df['Seq'].str[:-6]
    # identify non NNS sequences
    df['NNS'] = df['Seq'].apply(dna.is_bias_on)
    # identify seq with STOP codons
    df['STOP'] = df['Seq'].apply(dna.code_4_any)
    # filter by previous criteria and translate to AA
    dfc = df[(df['NNS']==True) & (df['STOP']==False)]
    # Error if any seq complain the filters
    assert(dfc.shape[0]>0)
    dfc['aa'] = dfc['Seq'].apply(dna.translate2aa)
    # save?
    # For each unique sequences group
    # sum sequences coding for the same
    grp_df = dfc.groupby(['aa'], as_index=False)['Reads'].agg(sum)
    # grp_df.rename(columns={'Seq': 'na'}, inplace=True)
    grp_df.rename(columns={'aa': 'Seq'}, inplace=True)

    if norm:
        diversity = grp_df['Seq'].unique().shape[0]
        grp_df['nReads'] = ( grp_df['Reads'] / grp_df['Reads'].sum() ) * diversity
        return grp_df[[ 'Seq', 'Reads', 'nReads']]
    else:
        return grp_df[['Seq', 'Reads']]

    # Data for Shannon entroy
    # transformed = list()
    # peptides = five4c['aa'].unique()

    # for pep in peptides:
        # coding = five4c[five4c['aa'] == pep]

        # coding['freq'] = coding['counts'] /  coding['counts'].sum()
        # variants = get_variants(pep)
        # suma = list()
        # for idx,c in coding.iterrows():
        #    suma.append(c['freq']*np.log2(c['freq']))
        # entropy = -1*(1.0/np.log2(variants))*sum(suma)
        # total_count = coding['Reads'].sum()
        # print pep, total_count, entropy
        # transformed.append([pep, total_count, entropy,coding.shape[0]])
        # transformed.append([pep, total_count, entropy,coding.shape[0]])

    # return pd.DataFrame(transformed,columns=['aa','counts','e','variants'])
    return grp_df


def b2h_random_variants(fasta_file, norm=True):
    """Pipeline for preprocess of random libraries from b2h.

        Read demultiplexed fasta files, filter non NNS and Frameshift sequences.
    Returns a Dataframe with the sequences group by AA

    Parameters
    ----------
        fasta_file : str
            Path to the fasta file

    Returns
    -------
        Dataframe
            Return a dataframe with filtered sequneces
            ['referenceId', 'Seq', 'Reads', 'nReads']


    """
    # load data
    df = read_fasta(fasta_file, norm=False, translate2aa=False)
    # remove extra stop codons
    # df['na'] = df['Seq'].str[:-6]
    # identify non NNS sequences
    df['NNS'] = df['Seq'].apply(dna.is_bias_on)
    # identify seq with STOP codons
    df['STOP'] = df['Seq'].apply(dna.code_4_any)
    # filter by previous criteria and translate to AA
    dfc = df[(df['NNS']==True) & (df['STOP']==False)]
    # Error if any seq complain the filters
    assert(dfc.shape[0]>0)
    dfc['aa'] = dfc['Seq'].apply(dna.translate2aa)
    # group and sum reads
    grup_sum = dfc.groupby(['aa'],as_index=False)['Reads'].agg(sum)
    # group and count variants
    grup_variants = dfc.groupby(['aa'],as_index=False)['Seq'].agg(len)
    # merge and change name
    d = pd.merge(grup_sum, grup_variants ,on='aa')
    d.rename(columns={'Seq':'Vars'},inplace=True)
    print(fasta_file)
    return d


def b2h_random_entropy(fasta_file, norm=True):
    """Pipeline for preprocess of random libraries from b2h.

        Read demultiplexed fasta files, filter non NNS and Frameshift sequences.
    Returns a Dataframe with the sequences group by AA

    Parameters
    ----------
        fasta_file : str
            Path to the fasta file

    Returns
    -------
        Dataframe
            Return a dataframe with filtered sequneces
            ['referenceId', 'Seq', 'Reads', 'nReads']


    """
    # load data
    df = read_fasta(fasta_file, norm=False, translate2aa=False)
    # remove extra stop codons
    # df['na'] = df['Seq'].str[:-6]
    # identify non NNS sequences
    df['NNS'] = df['Seq'].apply(dna.is_bias_on)
    # identify seq with STOP codons
    df['STOP'] = df['Seq'].apply(dna.code_4_any)
    # filter by previous criteria and translate to AA
    dfc = df[(df['NNS']==True) & (df['STOP']==False)]
    # Error if any seq complain the filters
    assert(dfc.shape[0]>0)
    dfc['aa'] = dfc['Seq'].apply(dna.translate2aa)
    # group and sum reads
    df_groups_aa =dfc.groupby('aa')
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
                parsed_aa.append([aa, df_grp['Reads'].sum(),df_grp.shape[0],  0])
        else:
            N = dna.possible_encondings(aa)
            df_grp['prob'] = df_grp['Reads'] / df_grp['Reads'].sum(skipna=True)
            df_grp['eprob'] = df_grp['prob'] * (np.log2( df_grp['prob']))

            shannon_entropy = -1 * ( (1.0/(np.log2(N))) *  (df_grp['eprob'].sum(skipna=True)))

            # Seq, Reads, Entropy
            parsed_aa.append([aa, df_grp['Reads'].sum(),df_grp.shape[0], shannon_entropy])
    return pd.DataFrame(parsed_aa, columns=['Seq', 'Reads','Var' ,'E'])


