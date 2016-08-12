#from __future__ import print_function
from dna_util import translate2aa
import argparse
import sys
import operator

__author__ = "C. corbi-Verge"
__email__ = "carles.corbi@gmail.com"
__about__ = ''' This is a small script to count and rank reads from a
             fasta file. This script is agnostic, do not compare or map
             the reads with any source or library'''

sample_id = sys.argv[1]
c1 = sys.argv[2]
pep = sys.argv[3]
c2 =sys.argv[4]
reading_quality_cutoff = 30

result = open(".%s.trim" % (sample_id), "w")
no_complaint = 0
complaint = 0
with open(sample_id, 'r') as forward:
    for forward_id in forward:

        header_id = forward_id
        header_id = str(forward_id.partition(' ')[0]).strip()
        sequence, line3, quality = [next(forward) for _ in range(3)]

        forward_id = forward_id.partition(' ')

        #Translate the Quality to numbers
        qual = [ord(c)-33 for c in quality.rstrip("\n")]

        if True:
            sequence = sequence[len(c1):len(pep)+len(c2)]
            #remove the quality of the barcode and the constant region
            qual = qual[len(c1):len(pep)+len(c2)]
            # lencutoff = int(self.peplength)





        #split the remain seq in diff frames in order to find the reverse constant region
        #Not always the seq has the 48bp, that why we use this method.

        dseq_q = sum(qual)/float(len(qual))
        if dseq_q >= reading_quality_cutoff:

            print >> result, ">"+header_id+"_F_"+str(len(sequence))+"_"+str(dseq_q)
            print >> result, sequence
            complaint +=1
        else:

            no_complaint +=1


