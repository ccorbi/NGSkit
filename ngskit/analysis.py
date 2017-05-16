"""
Collection of tools for sequence analysis
"""
import sklearn.metrics
import sklearn.cluster
from sklearn.cluster import AgglomerativeClustering

from sklearn.metrics import  silhouette_score
from Bio.SubsMat import MatrixInfo

import pandas as pd
from multiprocessing import Pool
from functools import partial

#
# - PSSM
#

def get_pfm(sequences):
    """Read list of sequence  and return frequency position Matrix.

    Parameters
    ----------
    sequences : array_like
       array of list containing the sequence in str format

    Returns
    -------
    dict
       Matrix on a dict format  where the key is the position, for each position the value is other
       dict with the counts  per amino acid i.e {'A': 44, 'C': 3, }


    """
    matrix = dict()
    for seq in sequences:
        for idx,p in enumerate(seq):
            # logger.debug("%i %s",idx,p)
            if idx not in matrix:
                matrix[idx] = {p:1}
            else:
                matrix[idx][p] = 1 + matrix[idx].get(p, 0)

    return matrix


def get_pssm(sequences, pseudocounts= 1.5):
    """Generate Position-Specific Scoring Matrix PSSM for sequences.

    Parameters
    ----------
    sequences : array_like
        array of list containing the sequence in str format

    Returns
    -------
    pandas.DataFrame
         Position-Specific Scoring Matrix , index is aminoacid, each column is the position



    """

    amino = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C",
    "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W", "-"]

    frequencies = get_aa_freq(sequences)


    N = len(sequences) + pseudocounts

    # get first element of the dict without know keys
    seq_len = len(sequences[0])
    matrix = list()

    for a in amino:
        # row = [a]
        row = list()
        for p in range(seq_len):
            row.append((frequencies[p].get(a, 0.0 ) + (pseudocounts/ 20 ) )/N)
        matrix.append(row)

    # convert to pandas.Dataframe
    m = pd.DataFrame(matrix, index=amino,columns=list(range(seq_len)))

    return m



def get_pwm(sequences):
    """Generate Position Weights Matrix

    Parameters
    ----------
    sequences : array_like
        array of list containing the sequence in str format

    Returns
    -------
    pandas.DataFrame
         Position Probability Matrix , index is aminoacid, each column is the position



    """

    amino = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C",
             "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]

    frequencies = get_pfm(sequences)
    # logger.debug(frequencies)

    N = len(sequences)
    # get first element of the dict without know keys
    seq_len = len(sequences[0])
    matrix = list()
    for a in amino:
        # row = [a]
        row = list()
        for p in range(seq_len):
            prob = frequencies[p].get(a, 0.0)/N
            row.append(np.log2(prob/.05))
        matrix.append(row)

    # convert to pandas.Dataframe
    m = pd.DataFrame(matrix, index=amino,columns=list(range(seq_len)))
    # logger.debug(m)
    return m


def binding_score(seq, pwm):
    """Score a sequence using a PWM.

    Parameters
    ----------
    seq: str

    pwm: Pandas.DataFrame

    Returns
    -------

    float

    """

    score = list()
    for pos, symbol in enumerate(seq):
        score.append(pwm.at[symbol, pos])

    return sum(score)

def dist_PWM(pwm1, pwm2):
    """Euclidina distance in the PWM.

    $ D_{(a,b)} = \frac{1}{W} \sum_{i=1}^{w} \frac{1}{2} \sum (a_{i,L} - b_{i,L})^2 $

    Parameters
    ----------

    Returns
    -------

    """
    assert pwm1.shape == pwm2.shape

    w = len(pwm1.columns)

    colum = list()

    for c in pwm1.columns:
        rows = list()
        for r in pwm1.index:
            rows.append((pwm1.at[r,c]-pwm2.at[r,c])**2)
        colum.append(sum(rows)*.5)

    return sum(colum)/float(w)



def entropy(sequences):
"""return a array with entropy for each postion.

Parameters
----------
ppm: pandas.DataFrame()

Returns
-------
array

"""

ppm = get_ppm(sequences)
entropy_vec = np.empty(ppm.shape[1])
# for position
for idx, a in enumerate(ppm.columns):
    e = np.empty(ppm.shape[0])
    # get entropy for all aminoacids
    for jdx, i in enumerate(ppm.index.get_values()):
        prob = ppm.at[i,a]
        e[jdx] = prob*(np.log(prob)/np.log(20))

    entropy_vec[idx] = -1 * np.nansum(e)

return entropy_vec

def mutual_information(sequences, i, j):
    """To identify correlated positions in peptide alignments, we used
    mutual information for all possible position pairs. Mutual information is computed as:.
(j)
    MI = sum 1-20 sum 1-20 p(i,j) log {P(i,j)/ P1i P2}

    where P(i,j) stands for the probability of having amino acid i at one position together
    with amino acid j at the other position in a peptide. P1(i) is the probability of having amino
    acid i at one position and P2(j) the probability of having amino acid j at the other position.
    MI=0 corresponds to independence of the two positions, whereas larger values indicate that
    knowing which amino acids are found at one position gives some information about which ones
    are expected at the other position. One limitation of mutual information is that non-zero
    values are often expected to be present by chance, even for randomly generated peptides.
    We therefore used the mutual information P-value as a filter (other statistical measures
    such as the Z-scores could also be used). All P-values have been computed by randomly
    shuffling the amino acids within the alignment columns. A threshold of 0.001 has been
    used to define correlated positions.

    Parameters
    ----------
    ppm : dict
    i: int
    j: int

    Returns
    -------
    float

    """
    aminoacids = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C",
    "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]
    # get Total seq
    #N = sum([X for X in ptm[i].values()])
    # n = sum([ sum(position.values()) for position in pfm.values()])
    frequencies = get_pfm(sequences)
    N = float(len(sequences))

    # transform to ....

    # insert frequency matrix, and positions1, position2
    # return MI of a pairf positions
    # given position K1, and poistion K2
    mi = list()
    for ai in aminoacids:
        for aj in aminoacids:
            mi.append(calc_mutual(sequences, frequencies, N, i, j ,  ai, aj ))
        #with Pool() as p:
        #    func = partial(calc_mutual,
        #               sequences=sequences,
        #               frequencies= frequencies,
        #               N = N,
        #               i=i,
        #               j=j,
        #               ai=ai)
        #    mi += p.map(func, aminoacids)


    return mi


def get_freq_pairs(sequences, i, j, ai, aj):

    filteri = [seq   for seq in sequences if seq[i] == ai]
    filterj =  [seq for seq in filteri if seq[j] == aj  ]

    return filterj

def calc_mutual(sequences, frequencies,N,  i, j ,  ai, aj ):

    prob_ai =  float(frequencies[i].get(ai, 0.0)) / N # get freq for that a in i
    prob_aj =  float(frequencies[j].get(aj, 0.0)) / N# get freq for that a in j
    freq_pairs_ij =  len(get_freq_pairs(sequences, i, j, ai, aj))
    prob_ij =  freq_pairs_ij / N# from fequency
    # print prob_ij, prob_ai, prob_aj, N, i, j , ai , aj
    try:

        r =  prob_ij*np.log(prob_ij/(prob_ai*prob_aj))

        #if r == np.float64('nan'):
        #    return 0.
        #else:
        return r

    except ZeroDivisionError:
        return np.nan

#
# CLUSTERING
#

def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]

def score_pairwise(seq1, seq2, matrix, gap_s, gap_e):
    score = 0
    gap = False
    for i in range(len(seq1)):
        pair = (seq1[i], seq2[i])
        if not gap:
            if '-' in pair:
                gap = True
                score += gap_s
            else:
                score += score_match(pair, matrix)
        else:
            if '-' not in pair:
                gap = False
                score += score_match(pair, matrix)
            else:
                score += gap_e
    return score

#seq1 = 'PAVKDLGAEG-ASDKGT--SHVVY----------TI-QLASTFE'
#seq2 = 'PAVEDLGATG-ANDKGT--LYNIYARNTEGHPRSTV-QLGSTFE'

#blosum = MatrixInfo.blosum62

#score_pairwise(seq1, seq2, blosum, -5, -1)

def find_clusters_(gpf, matrix ):

    m = matrix.as_matrix()
    big_score = 0.0
    ideal = 0
    for i in range(2,5):
        pred_kmean = sklearn.cluster.KMeans(n_clusters=i).fit_predict(m)
        # print 'gfp',gpf,'n',i,sklearn.metrics.silhouette_score(m,pred_kmean)
        cscore = sklearn.metrics.silhouette_score(m,pred_kmean)
        if cscore > big_score:
            ideal = i
            big_score = cscore



    return ideal

def generate_distance_matrix(dftemp, by='Seq', matrix=blosum):
    rowing = dict()
    for idx,row in dftemp.iterrows():

        coling = dict()

        for idj,r in dftemp.iterrows():

            #print r['aa'], score_pairwise(row['aa'],r['aa'],blosum,0,0)
            coling[r[by]] = score_pairwise(row[by],r[by],matrix,0,0)
        rowing[row[by]] = coling

    matrix = pd.DataFrame.from_dict(rowing)
    return matrix


def benchmark_clustering(matrix, limit=10, clustering_algorithm = AgglomerativeClustering, **Kargs ):
    """Return best Clustering sequences.

    Parameters
    ----------
    matrix = Pandas.DataFrame
        Pandas dataframe with all the information

    limit : int
        Maximum number of clusters

    clustering_algorithm : func
        by default skitlearn.Kmeans

    Returns
    -------

    labeled_matrix : pandas.DataFrame
        input matrix plus a Cluster label and dsitance to centroid  columns



    """
    #transform dataframe to matrix
    X = matrix.as_matrix()
    #init resutls
    cluster_labels = dict()
    # centroids = dict()
    silhouette_avg = dict()
    for i in range(2,limit):
        # tune perfomrmance using  n_jobs, or minibaches
        clusterer = clustering_algorithm(n_clusters=i, **Kargs)
        cluster_labels[i] = clusterer.fit_predict(X)
        # centroids[i] = clusterer.cluster_centers_
        # print 'gfp',gpf,'n',i,sklearn.metrics.silhouette_score(m,pred_kmean)
        silhouette_avg[i] = silhouette_score(X, cluster_labels[i])

    # get clustering with the lower score
    s = [k for k in sorted(silhouette_avg, key=silhouette_avg.get)]
    #transform values ToDo; figure out a much more elegant way to do it
    labeled_sequences = pd.DataFrame([matrix.index,cluster_labels[s[0]]]).T
    labeled_sequences.rename(columns={0:'Seq',1:'Cluster'}, inplace=True)
    labeled_matrix = pd.merge(matrix, labeled_sequences, left_index=True, right_on='Seq')

    # add distance to the centroid
    # labeled_matrix['Clust_distance'] = labeled_matrix.apply(distance_to_centroid, centroids=centroids[s[0]],
    #                                                  axis=1 )

    return labeled_matrix
