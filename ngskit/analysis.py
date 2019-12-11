# -*- coding: utf-8 -*-
"""
Collection of tools for sequence analysis
"""
import random
from multiprocessing import Pool
import math
from functools import partial
import subprocess
import pandas as pd
import numpy as np

import sklearn.metrics
import sklearn.cluster
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import  silhouette_score

from Bio.SubsMat import MatrixInfo

from  scipy.stats import entropy
from scipy.spatial.distance import cdist

from .utils.alphabets import *
from .utils.seqs import *


#
# - PSSM
#

def get_pfm(sequences):
    """Read list of sequence  and return counts position Matrix.

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


def get_ppm(sequences, pseudocounts= .0000001, elements=AMINOACIDS):
    """Generate Position-Specific Probability Matrix for sequences.

    Parameters
    ----------
    sequences : array_like
        array of list containing the sequence in str format

    Returns
    -------
    pandas.DataFrame
         Position-Specific Scoring Matrix , index is aminoacid, each column is the position



    """


    frequencies = get_pfm(sequences)


    N = len(sequences) + pseudocounts

    # get first element of the dict without know keys
    seq_len = len(sequences[0])
    matrix = list()

    for a in elements:
        # row = [a]
        row = list()
        for p in range(seq_len):
            row.append((frequencies[p].get(a, 0.0 ) + pseudocounts )/N)
        matrix.append(row)

    # convert to pandas.Dataframe
    m = pd.DataFrame(matrix, index=elements,columns=list(range(seq_len)))

    return m



def get_pwm(sequences,  pseudocounts= .0000001, elements=AMINOACIDS):
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



    frequencies = get_pfm(sequences)


    N = len(sequences) + pseudocounts

    background = 1/ len(elements)
    # get first element of the dict without know keys
    seq_len = len(sequences[0])
    matrix = list()

    for a in elements:
        # row = [a]
        row = list()
        for p in range(seq_len):
            prob = (frequencies[p].get(a, 0.0 ) + pseudocounts )/N
            row.append(np.log2(prob/background))
        matrix.append(row)

    # convert to pandas.Dataframe
    m = pd.DataFrame(matrix, index=elements,columns=list(range(seq_len)))

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
    """Euclidina distance between PWM.

    $ D_{(a,b)} = \frac{1}{W} \sum_{i=1}^{w} \frac{1}{2} \sum (a_{i,L} - b_{i,L})^2 $

    Parameters
    ----------

    Returns
    -------

    """
    # row = aa
    # col = position
    assert pwm1.shape == pwm2.shape

    w = len(pwm1.columns)

    colum = list()

    for c in pwm1.columns:
        rows = list()
        for r in pwm1.index:
            rows.append( (pwm1.at[r,c]-pwm2.at[r,c])**2) 
        colum.append(  np.sqrt(sum(rows)))

    return sum(colum)/float(w)


def ndist_PWM(pwm1, pwm2):
    """Normalized Euclidina distance introduced by Harbison et al 
    (Harbison et al 2004).  The score is between 1 complete dis-similar,
     and 0, perfect identity, and a distance of 0.20 indicates that 
     the two motifs differ by about 20%. We can consider to PWM to 
     represent the same motif is they score below 0.18
     
    """

    re = list()
    
    dist = cdist(pwm1.values.T,pwm2.values.T)
    for c in pwm1.columns:    
        re.append(dist[c][c])
    return (1/(np.sqrt(2)* 7))*sum(re)


def get_pem(sequences, alphabet, apply_correction=True):
    """return a matrix with IC for each postion and base/aminoacid.

    Parameters
    ----------
    sequences: array_like
        aligned sequences 
    alphabet: array_like 
        DNA, RNA or protein alphabet
    apply_correction: bool, optional
        the small-sample correction for an alignment of n letters

    Returns
    -------
    Pandas Dataframe
        Matrix contaning the the information content for each letter and position

    """

    ppm = get_ppm(sequences, pseudocounts=0)
    #entropy_vec = np.empty(ppm.shape[1])
    # small-sample correction
    n = len(sequences)
    if apply_correction:
        E = 1/np.log(2)
        ss_correction = E * ( (len(alphabet) - 1) /  (2*n ) ) 
    else:
        ss_correction = 0.0 


    # for position

    bits_df = freq2bits(ppm, ss_correction, alphabet)

    return bits_df

def freq2bits(df, ss_correction = 0, alphabet = ['A', 'C', 'G', 'T']):

    # each row is a position
    new_matrix  = list()
    for idx, c in enumerate(df.columns):
        h = entropy(df[c].values, base=2)
        # no correction here
        r = np.log2(len(alphabet)) - ( h + ss_correction )
        temp = df[c].values*r
        new_matrix.append(temp.tolist())
        #new_matrix = np.asarray(new_matrix)
    bits_df = pd.DataFrame(new_matrix, columns=alphabet)

    return bits_df.T


def mi_matrix(sequences, n_jobs=4):

    """

    """

    # get Total seq
    #N = sum([X for X in ptm[i].values()])
    # n = sum([ sum(position.values()) for position in pfm.values()])
    frequencies = get_pfm(sequences)
    N = len(sequences[0])
    # N = list(range(N))
    # transform to ....

    # insert frequency matrix, and positions1, position2
    # return MI of a pairf positions
    # given position K1, and poistion K2

    mi_matrix = pd.DataFrame(index=list(range(N)), columns=list(range(N)) )

    for i in range(N):
        for j in range(i,N):

            if i!=j:
                mutualinfo = mutual_information(sequences, frequencies, i, j )
                mi_matrix.at[i,j]= mutualinfo
                mi_matrix.at[j,i]= mutualinfo
            else:
                mi_matrix.at[i,j]= 0.0
                mi_matrix.at[j,i]= 0.0

        #with Pool() as p:
        #    func = partial(calc_mutual,
        #               sequences=sequences,
        #               frequencies= frequencies,
        #               N = N,
        #               i=i,
        #               j=j,
        #               ai=ai)
        #    mi += p.map(func, aminoacids)


    return mi_matrix


def get_freq_pairs(sequences, i, j, ai, aj):

    filteri = [seq  for seq in sequences if seq[i] == ai]
    filterj =  [seq for seq in filteri if seq[j] == aj  ]

    return filterj


def mutual_information(sequences, frequencies,  i, j  ):

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
    frequencies: pandas
    i: int
    j: int

    Returns
    -------
    float

    """




    N= float(len(sequences))

    suma = np.zeros([20,20])

    for idx,ax in enumerate(AMINOACIDS):
        for jdx,ay in enumerate(AMINOACIDS):
            # get prob aminoacid X on postion i
            prob_ai =  float(frequencies[i].get(ax, 0.0)) / N
            # get prob aminoacid Y on position j
            prob_aj =  float(frequencies[j].get(ay, 0.0)) / N
            # prob of co-appear
            freq_pairs_ij =  len(get_freq_pairs(sequences, i, j, ax, ay))
            prob_ij =  freq_pairs_ij / N# from fequency

            try:

                r =  prob_ij* ( np.log(prob_ij/(prob_ai*prob_aj)))


                suma[idx][jdx] = r

            except ZeroDivisionError:
                suma[idx][jdx] = 0.0

    return np.nansum(suma)


def mutual_informationm(sequences, frequencies,  i, j , n_jobs=8 ):

    """To identify correlated positions in peptide alignments, we used
    mutual information for all possible position pairs. Mutual information is computed as:.
(j)
    MI = sum 1-20 sum 1-20 p(i,j) log {P(i,j)/ P1i P2j}

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
    frequencies: pandas
    i: int
    j: int

    Returns
    -------
    float

    """


    N= float(len(sequences))

    suma = np.zeros([20,20])

    for idx,ax in enumerate(AMINOACIDS):
        for jdx,ay in enumerate(AMINOACIDS):
            # get prob aminoacid X on postion i
            prob_ai =  float(frequencies[i].get(ax, 0.0)) / N
            # get prob aminoacid Y on position j
            prob_aj =  float(frequencies[j].get(ay, 0.0)) / N
            # prob of co-appear
            freq_pairs_ij =  len(get_freq_pairs(sequences, i, j, ax, ay))
            prob_ij =  freq_pairs_ij / N# from fequency

            try:

                r =  prob_ij* ( np.log(prob_ij/(prob_ai*prob_aj)))


                suma[idx][jdx] = r

            except ZeroDivisionError:
                suma[idx][jdx] = 0.0

    return np.nansum(suma)

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

def generate_distance_matrix(dftemp, by='Seq', matrix=MatrixInfo.blosum62):
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



def pwm_mean(pwm):
    b = .05
    sx = .0
    for col in pwm.columns:
        for aa in pwm.index:
            logodds = pwm.at[aa, col]
            if math.isnan(logodds):
                continue
            if math.isinf(logodds) and logodds < 0:
                continue
            
            p = b * math.pow(2, logodds)
            sx += p * logodds
    return sx

def pwm_std(pwm):
    b =.05
    variance = 0.0
    for i in pwm.columns:
        sx = 0.0
        sxx = 0.0
        for letter in pwm.index:

            logodds = pwm.at[letter, i]
            if math.isnan(logodds):
                continue
            if math.isinf(logodds) and logodds < 0:
                continue

            p = b * math.pow(2, logodds)
            sx += p * logodds
            sxx += p * logodds * logodds
        sxx -= sx * sx
        variance += sxx
    variance = max(variance, 0)  # to avoid roundoff problems
    return math.sqrt(variance)


class ScoreDistribution(object):
    """Class representing approximate score distribution for a given motif.
    Utilizes a dynamic programming approach to calculate the distribution of
    scores with a predefined precision. Provides a number of methods for calculating
    thresholds for motif occurrences.
    """

    def __init__(self, pwm=None, precision=10 ** 3):
        """Initialize the class."""
        self.min_score = sum(pwm.min())
        self.interval = sum(pwm.max()) - self.min_score
        self.n_points = precision * len(pwm.columns)
        self.ic = pwm_mean(pwm)
        self.step = self.interval / (self.n_points - 1)
        self.mo_density = [0.0] * self.n_points
        self.mo_density[-self._index_diff(self.min_score)] = 1.0
        self.bg_density = [0.0] * self.n_points
        self.bg_density[-self._index_diff(self.min_score)] = 1.0

        
        # cal
        for pos in pwm.columns:
            mo_new = [0.0] * self.n_points
            bg_new = [0.0] * self.n_points
            for letter in pwm.index:
                bg = .05
                score = pwm.at[letter, pos]
                mo = pow(2, pwm.at[letter, pos]) * bg
                d = self._index_diff(score)
                for i in range(self.n_points):
                    mo_new[self._add(i, d)] += self.mo_density[i] * mo
                    bg_new[self._add(i, d)] += self.bg_density[i] * bg
            self.mo_density = mo_new
            self.bg_density = bg_new


    def _index_diff(self, x, y=0.0):
        #print(int((x - y + 0.5 * self.step) / self.step))
        return int((x - y + 0.5 * self.step) / self.step)

    def _add(self, i, j):
        return max(0, min(self.n_points - 1, i + j))

    def modify(self, scores, mo_probs, bg_probs):
        mo_new = [0.0] * self.n_points
        bg_new = [0.0] * self.n_points
        for k, v in scores.items():
            d = self._index_diff(v)
            for i in range(self.n_points):
                mo_new[self._add(i, d)] += self.mo_density[i] * mo_probs[k]
                bg_new[self._add(i, d)] += self.bg_density[i] * bg_probs[k]
        self.mo_density = mo_new
        self.bg_density = bg_new

    def threshold_fpr(self, fpr):
        """Approximate the log-odds threshold which makes the type I error (false positive rate)."""
        i = self.n_points
        prob = 0.0
        while prob < fpr:
            i -= 1
            prob += self.bg_density[i]
        return self.min_score + i * self.step

    def threshold_fnr(self, fnr):
        """Approximate the log-odds threshold which makes the type II error (false negative rate)."""
        i = -1
        prob = 0.0
        while prob < fnr:
            i += 1
            prob += self.mo_density[i]
        return self.min_score + i * self.step

    def threshold_balanced(self, rate_proportion=1.0, return_rate=False):
        """Approximate log-odds threshold making FNR equal to FPR times rate_proportion."""
        i = self.n_points
        fpr = 0.0
        fnr = 1.0
        while fpr * rate_proportion < fnr:
            i -= 1
            fpr += self.bg_density[i]
            fnr -= self.mo_density[i]
        if return_rate:
            return self.min_score + i * self.step, fpr
        else:
            return self.min_score + i * self.step

    def threshold_patser(self):
        """Threshold selection mimicking the behaviour of patser (Hertz, Stormo 1999) software.
        It selects such a threshold that the log(fpr)=-ic(M)
        note: the actual patser software uses natural logarithms instead of log_2, so the numbers
        are not directly comparable.
        """
        return self.threshold_fpr(fpr=2 ** -self.ic)




#            generate_logo(df['Seq'], filename='./temp/3_iterative_{}_{:.5f}_{}'.format(n,q,df.shape[0]))
def generate_logo(sequences, seq_len=80, filename='designs', **kwargs  ):
    """quick logo generation.

    Parameters
    ----------
    sequences : array_like

    seq_len : int


    filename : str

    Returns
    -------

    """
    #--label
    # annotate = xlabel
    #format == fformat
    #                 'title': '', 
    #             'annotate':'',
    # my defaults
    options = {'fineprint':'',
             'format':'PNG',
             'composition':'None' ,
             'units':'bits',
             'color-scheme':'chemistry',
             'stacks-per-line':300 }
    options.update(kwargs)
    
    ohandler = open(filename + '.fasta', 'w')
    for seq in sequences:
        print(">{}".format(seq), file=ohandler)
        print("{}".format(seq), file=ohandler)
    ohandler.close()
    
    base_commd = 'weblogo -f {} -o {}'.format(filename + '.fasta',
                                              filename + '.{}'.format(options['format']))
    # apply kwargs here
    for label, data in options.items():
        if data:
            base_commd = base_commd + f' --{label} {data} '

    
    
    command = subprocess.Popen(base_commd.split())
    if command.wait() != 0:
        output, error = command.communicate()
        print('Error Generating the logo')
        print(error, output)
        print(base_commd)

    return      



# this comes from Bio.motif, need reimplementation
# I am considering to rewrite the whole thing



# Default scales
# The following references were used for the default physicochemical scales:
# • Aliphatic index:
# Ikai AJ. Thermostability and aliphatic index of globular proteins. J Biochem 88, 1895-1898
# (1980).
# 5
# • Bulkiness scale:
# Zimmerman JM, Eliezer N, Simha R. The characterization of amino acid sequences in proteins
# by statistical methods. J Theor Biol 21, 170-201 (1968).
# • Hydrophobicity scale:
# Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein.
# J Mol Biol 157, 105-32 (1982).
# • pK values:
# http://emboss.sourceforge.net/apps/cvs/emboss/apps/iep.html
# • Polarity scale:
# Grantham R. Amin