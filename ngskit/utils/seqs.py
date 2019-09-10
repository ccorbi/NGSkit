import random
from  .alphabets import *


# Utils

def guess_alphabet(seq):

    coding_letters = set(seq)

    if len(coding_letters) > 4:

        return AMINOACIDS
    
    if is_aa_specific(coding_letters):
        return AMINOACIDS

    return BASES
         

def is_aa_specific(seq):

    AA = ["R", "H", "K", "D", "E", "S", "N", "Q",  "P", "V", "I", "L", "M", "F", "Y", "W", ]


    for i in seq:
        if i in AA:
            return True
    else:
        return False

# Generate Sequences


def rand_peptide(n=7):

    rpep = list()

    for i in range(n):
        rpep.append(random.choice(AMINOACIDS))

    return ''.join(rpep)


def rand_nucleotide(n=7):

    # need to add GC content control

    rpep = list()

    for i in range(n):
        rpep.append(random.choice(BASES))

    return ''.join(rpep)



# Chou-Fasman algorithm
# http://crdd.osdd.net/raghava/betatpred/chou.html
# http://crdd.osdd.net/raghava/betatpred/ref.html#ref12
# Conformational parameters and positional frequencies for helix,ß-sheet and ß-turn residues.
chou_fasman_data = [ ['A', 142,   83,   66,   0.06,   0.076,  0.035,  0.058]
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

chou_fasman_df = pd.DataFrame(chou_fasman_data, columns=['Aa','Sa','Sb', 'St', 'f0', 'f1','f2','f3' ])

def generate_peptide_choufasman(n,le=10, label = 'H'):
    
    probs = {'H':'Pa','E':'Pb','T':'Pt'}
    random_cf = list()
    for i in range(n):
        p = [chou_fasman_df.sample(weights=probs[label], axis=0)['Aa'].get_values()[0] for k in range(le)]
        random_cf.append(''.join(p))
        
    return random_cf


# Simple physical properties

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


# ToDO: move it to a JSON maybe more elegant

hydrophobicity: {"A": 1.80000,
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
                                "V": 4.20000,}
aa_volume: {"A":91.5000  ,
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
                            "V":141.7000 ,}
aa_polarity: {"A":0.0000, 
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
                            "V":0.1300, }
aa_pK:{"A":0.0000, 
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
                            "V":0.0000, }
aa_hydrophilicity:{
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
                            "V":  -1.5000,}

aa_asa:{"A" :27.8000, 
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
                            "V" :23.7000,}

aa_local_flexibility:{"A":705.42000, 
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
                            "V":4474.4199,}

aa_mass:{"A":70.079  ,
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
                            "V":99.133  ,    }

aa_charge: {"A": 0,
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
                                "V": 0,
                                
}

aa_Molweights = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
           'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
           'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
           'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }

#Hheisenberg
aa_heHydrophobicity = {
                'A':  0.620,  'R': -2.530  , 'N': -0.780  ,'D': -0.900  ,
                'C':  0.290  ,'Q': -0.850  ,'E': -0.740  ,'G':  0.480  ,
                'H': -0.400  ,'I':  1.380  , 'L':  1.060  ,  'K': -1.500  ,
                'M':  0.640  ,'F':  1.190  , 'P':  0.120  ,'S': -0.180  ,
                'T': -0.050  ,'W':  0.810  , 'Y':  0.260  , 'V':  1.080}

# kdHydrophobicity
#A simple method for displaying the hydropathic character of a protein. Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32.
aa_kdHydrophobicity = {'I':4.5	,
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

aa_ttHydrophobicity = {'D':  -3.27, 
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

## aggregated properties

aa_property_tables= {'hydrophobicity': aa_hydrophobicity,
                    'volume': aa_volume,
                    'polarity':aa_polarity, 
                    "pK":, aa_pK,
                    "hydrophilicity":,aa_hydrophilicity,
                    "asa":aa_asa,
                    "local_flexibility":aa_local_flexibility,
                    "mass":aa_mass,
                    "charge": aa_charge 

}