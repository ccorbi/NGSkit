import random
from  ngskit.utils.alphabets import *

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
