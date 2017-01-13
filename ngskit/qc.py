from __future__ import print_function
from __future__ import division

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
