# -*- coding: utf-8 -*-
"""
Collection of tools for sequence analysis
"""
from __future__ import print_function
from __future__ import division
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

AMINOACIDS = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C",
    "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W", ]


THREE_LETTERS = {"ARG":"R", "HIS":"H", 
                "LYS":"K", "ASP": "D",
                "GLU": "E","SER": "S",
                "THR": "T", "ASN": "N",
                "GLN" :"Q", "CYS" :"C",
                "GLY": "G", "PRO" : "P",
                "ALA": "A", "VAL": "V", 
                "ILE": "I", "LEU": "L", 
                "MET":  "M", "PHE": "F", 
                "TYR" :"Y", "TRP": "W" }

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



def entropy(sequences):
    """return a array with entropy for each postion.

    Parameters
    ----------
    ppm: pandas.DataFrame()

    Returns
    -------
    array

    """

    ppm = get_ppm(sequences, pseudocounts=0)
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





#
# 
#

def rand_peptide(n=7):

    rpep = list()

    for i in range(n):
        rpep.append(random.choice(AMINOACIDS))

    return ''.join(rpep)


data = [ ['A', 142,   83,   66,   0.06,   0.076,  0.035,  0.058]
,['R',  98,   93,   95,   0.070,  0.106,  0.099,  0.085]
,['N', 101,   54,  146,   0.147,  0.110,  0.179,  0.081]
, ['D',  67,   89,  156,   0.161,  0.083,  0.191,  0.091]
, ['C',  70,  119,  119,   0.149,  0.050,  0.117,  0.128]
, ['E', 151,   37,   74,   0.056,  0.060,  0.077,  0.064]
,['Q', 111,  110,   98,   0.074,  0.098,  0.037,  0.098]
, ['G',  57,   75,  156,   0.102,  0.085,  0.190,  0.152]
, ['H', 100,   87,   95,   0.140,  0.047,  0.093,  0.054]
, ['I', 108,  160,   47,   0.043,  0.034,  0.013,  0.056]
, ['L', 121,  130,   59,   0.061,  0.025,  0.036,  0.070]
, ['K', 114,   74,  101,   0.055,  0.115,  0.072,  0.095]
, ['M', 145,  105,   60,   0.068,  0.082,  0.014,  0.055]
, ['F', 113,  138,   60,   0.059,  0.041,  0.065,  0.065]
, ['P',  57,   55,  152,   0.102,  0.301,  0.034,  0.068]
, ['S',  77,   75,  143,   0.120,  0.139,  0.125,  0.106]
,['T',  83,  119,   96,   0.086,  0.108,  0.065,  0.079]
,['W', 108,  137,   96,   0.077,  0.013,  0.064,  0.167]
, ['Y',  69,  147,  114,   0.082,  0.065,  0.114,  0.125]
,['V', 106,  170,   50,   0.062,  0.048,  0.028,  0.053]]

choufasman = pd.DataFrame(data, columns=['Aa','Sa','Sb', 'St', 'f0', 'f1','f2','f3' ])


def generate_choufus(n,le=10, label = 'H'):
    
    probs = {'H':'Pa','E':'Pb','T':'Pt'}
    random_cf = list()
    for i in range(n):
        p = [choufasman.sample(weights=probs[label], axis=0)['Aa'].get_values()[0] for k in range(le)]
        random_cf.append(''.join(p))
        
    return random_cf


#
# get simple physical properties

def net_charge(seq):
# DEPRICATED    
    POS = ['R','K']
    NEG = ['E','D']
    charge = 0
    for s in seq:
        if s in POS:
            charge +=1
        if s in NEG:
            charge -=1
    return charge

# move 
aa_property_tables= {'hydrophobicity': {
                                "A": 1.80000,
                                "R": -4.5000,
                                "N": -3.5000,
                                "D": -3.5000,
                                "C": 2.50000,
                                "Q": -3.5000,
                                "E": -3.5000,
                                "G": -0.4000,
                                "H": -3.2000,
                                "I": 4.50000,
                                "L": 3.80000,
                                "K": -3.9000,
                                "M": 1.90000,
                                "F": 2.80000,
                                "P": -1.6000,
                                "S": -0.8000,
                                "T": -0.7000,
                                "W": -0.9000,
                                "Y": -1.3000,
                                "V": 4.20000,},
                    'volume': {
                            "A":91.5000  ,
                            "R":202.0000 ,
                            "N":135.2000 ,
                            "D":124.5000 ,
                            "C":118.0000 ,
                            "Q":161.1000 ,
                            "E":155.1000 ,
                            "G":66.40000 ,
                            "H":167.3000 ,
                            "I":168.8000 ,
                            "L":167.9000 ,
                            "K":171.3000 ,
                            "M":170.8000 ,
                            "F":203.4000 ,
                            "P":129.3000 ,
                            "S":99.10000 ,
                            "T":122.1000 ,
                            "W":237.6000 ,
                            "Y":203.6000 ,
                            "V":141.7000 ,},
                    'polarity': {
                            "A":0.0000, 
                            "R":52.000, 
                            "N":3.3800, 
                            "D":40.700, 
                            "C":1.4800, 
                            "Q":3.5300, 
                            "E":49.910, 
                            "G":0.0000, 
                            "H":51.600, 
                            "I":0.1500, 
                            "L":0.4500, 
                            "K":49.500, 
                            "M":1.4300, 
                            "F":0.3500, 
                            "P":1.5800, 
                            "S":1.6700, 
                            "T":1.6600, 
                            "W":2.1000, 
                            "Y":1.6100, 
                            "V":0.1300, },
                    "pK":{
                            "A":0.0000, 
                            "R":12.480, 
                            "N":0.0000, 
                            "D":3.6500, 
                            "C":8.1800, 
                            "Q":0.0000, 
                            "E":4.2500, 
                            "G":0.0000, 
                            "H":6.0000, 
                            "I":0.0000, 
                            "L":0.0000, 
                            "K":10.530, 
                            "M":0.0000, 
                            "F":0.0000, 
                            "P":0.0000, 
                            "S":0.0000, 
                            "T":0.0000, 
                            "W":0.0000, 
                            "Y":10.700, 
                            "V":0.0000, },
                    "hydrophilicity":{
                            "A":  -0.5000, 
                            "R":  3.00000, 
                            "N":  0.20000, 
                            "D":  3.00000, 
                            "C":  -1.0000, 
                            "Q":  0.20000, 
                            "E":  3.00000, 
                            "G":  0.00000, 
                            "H":  -0.5000, 
                            "I":  -1.8000, 
                            "L":  -1.8000, 
                            "K":  3.00000, 
                            "M":  -1.3000, 
                            "F":  -2.5000, 
                            "P":  0.00000, 
                            "S":  0.30000, 
                            "T":  -0.4000, 
                            "W":  -3.4000, 
                            "Y":  -2.3000, 
                            "V":  -1.5000,},
                    "asa":{   
                            "A" :27.8000, 
                            "R" :94.7000, 
                            "N" :60.1000, 
                            "D" :60.6000, 
                            "C" :15.5000, 
                            "Q" :68.7000, 
                            "E" :68.2000, 
                            "G" :24.5000, 
                            "H" :50.7000, 
                            "I" :22.8000, 
                            "L" :27.6000, 
                            "K" :103.000, 
                            "M" :33.5000, 
                            "F" :25.5000, 
                            "P" :51.5000, 
                            "S" :42.0000, 
                            "T" :45.0000, 
                            "W" :34.7000, 
                            "Y" :55.2000, 
                            "V" :23.7000,},
                    "local_flexibility":{
                            "A":705.42000, 
                            "R":1484.2800, 
                            "N":513.46010, 
                            "D":34.960000, 
                            "C":2412.5601, 
                            "Q":1087.8300, 
                            "E":1158.6600, 
                            "G":33.180000, 
                            "H":1637.1300, 
                            "I":5979.3701, 
                            "L":4985.7300, 
                            "K":699.69000, 
                            "M":4491.6602, 
                            "F":5203.8599, 
                            "P":431.96000, 
                            "S":174.76000, 
                            "T":601.88000, 
                            "W":6374.0698, 
                            "Y":4291.1001, 
                            "V":4474.4199,},
                    "mass":{
                            "A":70.079  ,
                            "R":156.188 ,
                            "N":114.104 ,
                            "D":115.089 ,
                            "C":103.144 ,
                            "Q":128.131 ,
                            "E":129.116 ,
                            "G":57.052  ,
                            "H":137.142 ,
                            "I":113.160 ,
                            "L":113.160 ,
                            "K":128.174 ,
                            "M":131.198 ,
                            "F":147.177 ,
                            "P":97.177  ,
                            "S":87.078  ,
                            "T":101.105 ,
                            "W":186.213 ,
                            "Y":163.170 ,
                            "V":99.133  ,    },
                    "charge": {
                                "A": 0,
                                "R": 1,
                                "N": 0,
                                "D": -1,
                                "C": 0,
                                "Q": 0,
                                "E": -1,
                                "G": 0,
                                "H": 0,
                                "I": 0,
                                "L": 0,
                                "K": 1,
                                "M": 0,
                                "F": 0,
                                "P": 0,
                                "S": 0,
                                "T": 0,
                                "W": 0,
                                "Y": 0,
                                "V": 0,},        }








# ToDO: move it to a JSON maybe more elegant

Molweights = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
           'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
           'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
           'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }

#Hheisenberg
heHydrophobicity = {
                'A':  0.620,  'R': -2.530  , 'N': -0.780  ,'D': -0.900  ,
                'C':  0.290  ,'Q': -0.850  ,'E': -0.740  ,'G':  0.480  ,
                'H': -0.400  ,'I':  1.380  , 'L':  1.060  ,  'K': -1.500  ,
                'M':  0.640  ,'F':  1.190  , 'P':  0.120  ,'S': -0.180  ,
                'T': -0.050  ,'W':  0.810  , 'Y':  0.260  , 'V':  1.080}

# kdHydrophobicity
#A simple method for displaying the hydropathic character of a protein. Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32.
kdHydrophobicity = {'I':4.5	,
                    'V':4.2	,
                    'L':3.8	,
                    'F':2.8	,
                    'C':2.5	,
                    'M':1.9	,
                    'A':1.8	,
                    'G':-0.4,	
                    'T':-0.7,	
                    'S':-0.8,	
                    'W':-0.9,	
                    'Y':-1.3,	
                    'P':-1.6,	
                    'H':-3.2,	
                    'E':-3.5,	
                    'Q':-3.5,	
                    'D':-3.5,	
                    'N':-3.5,	
                    'K':-3.9,	
                    'R':-4.5,}


# wwHydrophobicity
#Experimentally determined hydrophobicity scale for proteins at membrane interfaces. Wimley WC, White SH. Nat Struct Biol. 1996 Oct;3(10):842-8. Attribute assignment file .txt.
# https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/defineattrib/wwHydrophobicity.txt

# hhHydrophobicity
#Recognition of transmembrane helices by the endoplasmic reticulum translocon. Hessa T, Kim H, Bihlmaier K, Lundin C, Boekel J, Andersson H, Nilsson I, White SH, von Heijne G. Nature. 2005 Jan 27;433(7024):377-81, supplementary data. Attribute assignment file hhHydrophobicity.txt. In this scale, negative values indicate greater hydrophobicity.
 # https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/defineattrib/hhHydrophobicity.txt

# mfHydrophobicity
#Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers. Moon CP, Fleming KG. Proc Natl Acad Sci USA. 2011 Jun 21;108(25):10174-7, supplementary data. Attribute assignment file mfHydrophobicity.txt. In this scale, negative values indicate greater hydrophobicity.
# https://www.cgl.ucsf.edu/chimera/docs/ContributedSoftware/defineattrib/mfHydrophobicity.txt

# ttHydrophobicity
#An amino acid “transmembrane tendency” scale that approaches the theoretical limit to accuracy for prediction of transmembrane helices: relationship to biological hydrophobicity. Zhao G, London E. Protein Sci. 2006 Aug;15(8):1987-2001. Attribute assignment file ttHydrophobicity.txt (contributed by Shyam M. Saladi).

ttHydrophobicity = {'D':  -3.27, 
                    'E':  -2.90, 
                    'N':  -1.62, 
                    'Q':  -1.84, 
                    'K':  -3.46, 
                    'R':  -2.57, 
                    'H':  -1.44, 
                    'H':  -0.19, 
                    'P':  -1.44, 
                    'S':  -0.53, 
                    'T':  -0.32, 
                    'C':  -0.30, 
                    'M':  1.40,
                    'A':  0.38,
                    'V':  1.46,
                    'I':  1.97,
                    'L':  1.82,
                    'F':  1.98,
                    'W':  1.53,
                    'Y':  0.49}





#            generate_logo(df['Seq'], filename='./temp/3_iterative_{}_{:.5f}_{}'.format(n,q,df.shape[0]))
def generate_logo(sequences, seq_len=80, filename='designs', title = None, fineprint='',xlabel=None, fformat='PNG' ):
    """quick logo generation.

    Parameters
    ----------
    sequences : array_like

    seq_len : int


    filename : str

    Returns
    -------

    """
    # if pass , Folder name
    
    #--fineprint
    #--title
    #--label
    extras = {'fineprint':fineprint,
             'title': title, 
             'annotate':xlabel,
             'format':fformat}

    ohandler = open(filename + '.fasta', 'w')
    for seq in sequences:
        print(">{}".format(seq), file=ohandler)
        print("{}".format(seq), file=ohandler)
    ohandler.close()
    
    base_commd = 'weblogo -f {} -c chemistry -o {} -n {} -U bits --composition equiprobable '.format(filename + '.fasta',
                                                                            filename + '.{}'.format(fformat),
                                                                            seq_len)

    for label, data in extras.items():
        if data:
            base_commd = base_commd + f'--{label} {data} '

    
    
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