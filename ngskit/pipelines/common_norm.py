# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division
import sys
import glob
import pandas as pd
import numpy as np
import scipy as sp
# My package
from  ngskit.utils import dna


# Normalize Read counts

def norm_TC(df, scalar_factor = 1):
    """Very basic normalization. Total Counts. scalar_factor =  1

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with the raw reads.
    scalar_factor : int, float
        Integer or float, can be diversity of the smaples or the avg of read across
        the whole dataset

    Returns
    -------
    pandas.DataFrame
        with the colun nReads, Normilized reads

    """
    df['nReads'] = ( df['Reads'] / df['Reads'].sum() ) * scalar_factor

    return df

def norm_TC_AVG(df):
    """Very basic normalization. Total Counts.  scalar_factor= Avg of reads

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with the raw reads.
    scalar_factor : int, float
        Integer or float, can be diversity of the smaples or the avg of read across
        the whole dataset

    Returns
    -------
    pandas.DataFrame
        with the colun nReads, Normilized reads

    """
    scalar_factor = df['Reads'].mean(skipna=True)

    df['nReads'] = ( df['Reads'] / df['Reads'].sum() ) * scalar_factor

    return df


def norm_TC_COUNTS(df):
    """Very basic normalization. Total Counts. scalar_factor is the total of mapped peptides

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with the raw reads.
    scalar_factor : int, float
        Integer or float, can be diversity of the smaples or the avg of read across
        the whole dataset

    Returns
    -------
    pandas.DataFrame
        with the colun nReads, Normilized reads

    """
    scalar_factor = df['Seq'].unique().shape[0]
    df['nReads'] = ( df['Reads'] / df['Reads'].sum() ) * scalar_factor

    return df


def norm_deseq(df, read_columns):
    """Scale a dataframe using the deseq scaling.

    Uses :func:`scale_factor_deseq`

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with the raw reads.
    read_columns : list
        list with the columns with the read counts information

    Returns
    -------
    pandas.DataFrame
        with the colun nReads, Normilized reads
    """

    # rename columns to label norm reads
    nreads_colums = ['n'+col_name for col_name in read_columns]
    # get scalar factor to norm reads
    scale_factors = scale_factor_deseq(df, read_columns)
    # norm
    df[nreads_colums] = df[read_columns] / scale_factors

    return df


def scale_factor_deseq(df, read_columns):
    """Returns the scale factor according to he deseq paper. The columns of the
    dataframe are the samples.

    size factor :math:`\\hat{s}_{j}` for sample *j* (from DESeq paper).

    .. math::

        \\hat{s}_{j} = median_{i} (
        \\frac
            {k_{ij}}
            {
                \\left (
                \\prod_{v=1}^{m}
                    k_{iv}
                \\right )^{1/m}
           }
        )

    """
    # set nan values, if any, to zero
    df[read_columns].fillna('0', inplace=True)
    # calc the genes geometric mean over all samples
    gmean = df[read_columns].apply(sp.stats.gmean, axis=1)
    # keeps only the genes whose geometric mean is > 0
    gmean = gmean[gmean > 0]

    sample_factors = {}

    # calc the scaling factor for each sample
    for sample, seq in df[read_columns].iteritems():

        scale_factor = np.median(seq.loc[gmean.index] / gmean)

        sample_factors[sample] = scale_factor

    return pd.Series(sample_factors)

# I/O

def read_preprocessed(file_dir, norm=False, **Kargs):
    """read file with sequences already counted.

    Parameters
    ----------

    Returns
    -------
    Pandas DataFrame ['Seq','Reads'] no Normalized
    Pandas DataFrame ['Seq','Reads', 'nReads']  Normalized
    """

    df =  pd.read_csv(file_dir, names=['Seq','Reads'], **Kargs)
    # very basic normalization
    if norm:

        df = norm(df)

    return df

def read_fasta(fasta_file, norm=False, translate2aa=True, **Kargs):
    """Read a demultiplexed fasta file and return Dataframe.

    Parameters
    ----------

    Returns
    -------
    Pandas DataFrame ['Seq','Reads'] no Normalized
    Pandas DataFrame ['Seq','Reads', 'nReads']  Normalized
    """

    seqs = dict()
    with open(fasta_file, 'r') as input_file:
        for line in input_file:
            if line.startswith('>'):
                pass
            else:
                s = line.strip()

                if translate2aa:
                    aa = dna.translate2aa(s)

                    seqs[aa] = seqs.get(aa, 0) + 1
                else:

                    seqs[s] = seqs.get(s, 0) + 1

    df = pd.DataFrame.from_dict(data=seqs, orient='index')
    df.reset_index(inplace=True)
    df.rename(columns={'index':'Seq', 0:'Reads'}, inplace=True)

    if norm:
        # Apply normalization function

        df = norm(df)

    else:
        pass

    return df

def read_bowtie(filename, norm=False, **Kargs):
    """Read bowtie output.

    Parameters
    ----------

    Returns
    -------
    Pandas DataFrame ['Seq','Reads'] no Normalized
    Pandas DataFrame ['Seq','Reads', 'nReads']  Normalized


    """
    df = load_bowtieout(filename)
    # Join by seq, count and sort
    q = df.groupby(['referenceId'], as_index=False)['offset'].agg(len)
    q.sort_values(by='offset', ascending=False, inplace=True)
    q.rename(columns={'offset': 'Reads'}, inplace=True)

    if norm:
        # Apply  normalization function

        q = norm(q)

        return q[['referenceId', 'Reads','nReads']]
    else:
        return q[['referenceId', 'Reads']]


def read_blast(filename, aligment_length, norm=False, **Kargs):
    """Read blasted file.

    Parameters
    ----------
    filename : str
        Fasta output file dir

    aligment_length :

    norm : Function
        Normalization fucntion
    Returns
    -------
    Pandas DataFrame ['referenceId', 'queryId', 'Reads'] no Normalized
    Pandas DataFrame [''referenceId', 'queryId', 'Reads', 'nReads']  Normalized


    """
    df = load_blastout(filename)
    # filter cutoff identity, aligment length and mismatchCount
    mismatchCount = Kargs.get('mismatchCount', d=0)
    identity = Kargs.get('identity', d=100.0)

    df = df[(df['percIdentity'] >= identity)
            & (df['alnLength'] == aligment_length)
            & (df['mismatchCount'] <= mismatchCount)]

    # Join by seq, count and sort
    q = df.groupby(['referenceId'], as_index=False)['eVal'].agg(len)
    q.sort_values(by='eVal', ascending=False, inplace=True)
    q.rename(columns={'eVal': 'Reads'}, inplace=True)

    if norm:
        # Apply  normalization function
        # logg with type is appling
        q = norm(q)
        # q['nReads'] = ( q['Reads'] / q['Reads'].sum() ) * diversity

        return q[['referenceId', 'Reads','nReads']]
    else:
        return q[['referenceId',  'Reads']]

# I/O More generic

def load_blastout(filename):
    """load blast output to Dataframe.


    Parameters
    ----------

    Returns
    -------

    Notes
    -----
    .: How to generate output blast:

    blastall -p blastn -d design_lib.fasta -i smaples_id_demultiplex.fasta -o blastoutput -v1 -b1 -m8
        or
    blastn -db /home/kimlab2/ccorbi/optim_lib/Kim/fasta_libs/optm_library.fasta -query $i  -evalue 1 -out $i.out  -outfmt '6 qacc sacc pident length  mismatch gapopen qstart qend sstart send evalue bitscore'


    """
    df = pd.read_csv(filename, delim_whitespace=True, names=[
        'queryId',
        'referenceId',
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
    """load bowtie output to DataFrame.

    Parameters
    ----------
    filename : str
        bowtie output file dir




    Returns
    -------
    DataFrame
        Pandas dataframe with  ['queryId', 'strand', 'referenceId',
                                'off',  'Seq', 'otheInstances',
                                'quality', 'mismatchDescriptor' ]
    Notes
    -----
    .: Bowite for 0 mismathces output must be generated (bowtie 2 output) as
        bowtie --best -v X [library/path] -f [fasta/path] > [output/path]

    .: Bowtie library:
        bowtie-build -f [fasta/library/path] [name library]

    .: Historical
        bowtie -n 0 [library/path] -f [fasta/path] > [output/path]

    """
    df = pd.read_csv(filename, names=[
        'queryId',
        'strand',
        'referenceId',
        'offset',
        'Seq',
        'Phred',
        'otherInstances',
        'mismatchDescriptor',
        ], sep='\t')


    return df




#
# CUSTOM
# Calculations
#
def dropOutScore(row, sampling_points=[0,3,7,14], dbl_times = [0,2.4,4.5,10], normed=True, itemization=False, field='nReads_T{}'):
    """This fucntion calcs the drop off score for each peptide_id

    Parameters
    ----------
    :param Serie row: Row from a dataframe, or dict
    :param list sampling_points: Time points, [0,3,7,14]
    :param list dbl_times: Double time of the cell line. [0,2.4,4.5,10]
    :param bool normed: Divison all the time points by the biggest value in the time serie
    :param bool itemization: return slope for each time point and the dropscore
    :param int scalar: scalar to multiply dropscore [increse the signal]
    :param str field: reads, of Cell viability base

    Returns
    -------
        :type float: dropout score
    if itemization is True:
        :type list: Slopes for each time point
        :type float: dropout score

    """

    # Init variables
    sampling_points.sort()
    data = []
    # Get initial point
    T0 = row[field.format(0)]

    # sync time points and reads
    data_points = [row[field.format(t)] for t in sampling_points]

    # Divison all the time points by the biggest value in the time serie
    # Avoid bias by highread peptides, relative comparision

    if normed:
        Norm_vector = max(data_points)
    else:
        Norm_vector =T0
    # get number of points

    n = len(sampling_points)

    for i, t in enumerate(sampling_points):
        # Skip T = 0
        if i == 0:
            data.append(0)
        else:
            s = slope_i(y0= T0/Norm_vector,
                        y1 = row[field.format(t)]/Norm_vector,
                        deltaX = dbl_times[i])
            data.append(s)

    # Sum Slopes *  1 / N-1
    dropout = (sum(data)*(1/(float(n-1))))

    # return slope for each time point and the dropscore
    if itemization:
        return data,dropout*100
    else:
        return dropout*100

def dropOutScoreLog2(row, sampling_points=[0,3,7,14], dbl_times = [0,2.4,4.5,10], normed=True, itemization=False, field='nReads_T{}'):
    """This fucntion calcs the drop off score for each peptide_id

    Parameters
    ----------
    :param Serie row: Row from a dataframe, or dict
    :param list sampling_points: Time points, [0,3,7,14]
    :param list dbl_times: Double time of the cell line. [0,2.4,4.5,10]
    :param bool normed: Divison all the time points by the biggest value in the time serie
    :param bool itemization: return slope for each time point and the dropscore
    :param int scalar: scalar to multiply dropscore [increse the signal]
    :param str field: reads, of Cell viability base

    Returns
    -------
        :type float: dropout score
    if itemization is True:
        :type list: Slopes for each time point
        :type float: dropout score

    """

    # Init variables
    sampling_points.sort()
    data = []
    # Get initial point
    T0 = row[field.format(0)]

    # sync time points and reads
    data_points = [row[field.format(t)] for t in sampling_points]

    # Divison all the time points by the biggest value in the time serie
    # Avoid bias by highread peptides, relative comparision

    if normed:
        Norm_vector = max(data_points)
    else:
        Norm_vector =T0
    # get number of points

    n = len(sampling_points)

    for i, t in enumerate(sampling_points):
        # Skip T = 0
        if i == 0:
            data.append(0)
        else:
            s = slope_i(y0= np.log2(T0+.5),
                        y1 = np.log2(row[field.format(t)]+.5),
                        deltaX = dbl_times[i])
            data.append(s)

    # Sum Slopes *  1 / N-1
    dropout = (sum(data)*(1/(float(n-1))))

    # return slope for each time point and the dropscore
    if itemization:
        return data,dropout*100
    else:
        return dropout*100


def DP_linreg(row, sampling_points=[0,3,7,14], dbl_times = [0,2.4,4.5,10],field='nReads_T{}'):
    xis = list()
    yis = list()
    #m = row[['nReads_T0_RWP1','nReads_T4_RWP1','nReads_T7_RWP1','nReads_T14_RWP1']].max()
    for idx,t in enumerate(sampling_points):
        #for r in ['A', 'B', 'C']:
        xis.append(dbl_times[idx])
        yis.append(np.log2(row[field.format(t)]+.5))
    slope, i, r, p, e = sp.stats.linregress(xis,yis)

    return slope*100


def relative_increment(row, ref_value, target_value):
    """Relative increment.

    Parameters
    ----------

    Returns
    -------

    """
    increment = slope_i(y0 = row[ref_value],
                        y1 = row[target_value],
                        deltaX= row[target_value] + row[ref_value])


    return increment

def slope(y0, y1, x0, x1 ):
    """Resturn slope.

    Parameters
    ----------

    Returns
    -------

    """
    deltaX = float(x1) - float(x0)
    deltaY = float(y1) - float(y0)

    return deltaY/deltaX


def slope_i(y0, y1, deltaX):
    """Resturn slope.

    Parameters
    ----------

    Returns
    -------

    """
    deltaY = float(y1) - float(y0)

    return deltaY/float(deltaX)
