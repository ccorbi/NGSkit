from __future__ import print_function
from __future__ import division

import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


logger = logging.getLogger(__name__)


def gini(reads):
    """Gini index for a set of reads.

    Gini index, a common measure of income inequality in economics, can measure the evenness of sgRNA
    read counts [14]. It is perfectly normal for later time points in positive selection experiments
    to have higher Gini index since a few surviving clones (a few sgRNA with extreme high counts)
    could dominate the final pool while most of the other cells die (more sgRNAs with zero-count).
    In contrast, high Gini index in plasmid library, in early time points, or in negative selection
    experiments may indicate CRISPR oligonucleotide synthesis unevenness, low viral transfection
    efficiency, and over selection, respectively. At time 0 it should be between 0.1-0.2

    Parameters
    ----------
    reads : array_like

    Returns
    -------

    float

    """
    sorted_list = sorted(reads)
    height, area = 0, 0
    for value in sorted_list:
        height += value
        area += height - value / 2.
    fair_area = height * len(reads) / 2.
    return (fair_area - area) / fair_area


def read_fastq(inputfile):

    quality = list()
    with open(inputfile, 'r') as read1:

        for read1_id in read1:
            # Read 4 by 4
            # ID lane info, seq info etc
            # Read seq and Quality info

            read1_seq, read1_strand, read1_qual = [next(read1) for _ in range(3)]
            qual = [ord(c)-33 for c in read1_qual.strip()]
            quality.append(qual)

    return quality


class Stats(object):
    """Init dataframe to Save frequencies.

    Parameters
    ----------
    barcodes : list, iterable
        contains barcode objects to demultiplex sequneces

    Returns
    -------
    DataFrame

    idx       b1  c1  c2  b2 target_len1 ... target_lenX
    sample1
    sample2
    ...
    sampleX

    """
    def __init__(self, barcodes_list):
        super(Stats, self).__init__()
        
        index_ids = []
        columns_id = ['b1', 'c1', 'c2', 'b2']
        for barcode in barcodes_list:

            index_ids.append(barcode.id)
            if not str(barcode.trgt_len) in columns_id:
                columns_id.append(str(barcode.trgt_len))

        # init stats

        self._stats = pd.DataFrame(index=index_ids, columns=columns_id)
        self._stats.fillna(0, inplace=True)


    def write_freqs(self, read, barcode):
        """Update frequencies in the DataFrame.

        Parameters
        ----------
        read : dict
            Dict from the identification fuction with information about the match
        barcode : object

        Returns
        -------
        DataFrame

        """
        if read['b1'] != '':
            # Seq match barcode 1 (+1)
            self._stats.at[barcode.id, 'b1'] = self._stats.at[barcode.id,
                                                        'b1'] + 1
        if read['c1'] != '':
            # Seq match conc 1 (+1)
            self._stats.at[barcode.id, 'c1'] = self._stats.at[barcode.id,
                                                        'c1'] + 1
        if read['c2'] != '':
            # Seq match conc 2 (+1)
            self._stats.at[barcode.id, 'c2'] = self._stats.at[barcode.id,
                                                        'c2'] + 1
        if read['b2'] != '':
            # Seq match barcode 2 (+1)
            self._stats.at[barcode.id, 'b2'] = self._stats.at[barcode.id,
                                                        'b2'] + 1
        # Target sequnces, length +1
        if read['map']:
            if str(read['target_len']) in self._stats.columns:
                self._stats.at[barcode.id,
                            str(read['target_len'])] = self._stats.at[barcode.id,
                                                                    str(read['target_len'])] + 1

            else:
                self._stats.at[barcode.id, str(read['target_len'])] = 1
                self._stats.fillna(0, inplace=True)

    def save(self, name):
        
        self._stats.to_csv(name+'.csv')



def plot_readslens(stats_df,  sample_name):
   
    fig, ax = plt.subplots(figsize=(20,10))

    df = stats_df[stats_df['Sample']==sample_name]
    del df['Sample']
    df.columns = df.columns.astype(int,copy=False)
    df.sort_index(axis=1, inplace=True)
    ax.bar(df.columns,df.get_values()[0])
    ax.set_title(sample_name)
    plt.xticks(df.columns , df.columns, rotation='vertical')

    fig.savefig(sample_name+'_readlens.png')
    plt.close(fig)