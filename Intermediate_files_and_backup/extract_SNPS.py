import sys
from Bio import AlignIO
import pandas as pd
import numpy as np
import os
import operator



from Bio.Align import MultipleSeqAlignment


####READ ALN

aln=AlignIO.read(sys.argv[1].strip(),"fasta")



##select clade

aln2=MultipleSeqAlignment([])

for i in aln:
    if "U.3" in str(i.id):
        aln2.append(i)
        
##trim alignment 5'

for pos in range(aln2.get_alignment_length()):
    if not "-" in aln[:, pos]:
        p5p = pos
        break
        

##trim alignment 3'

for pos in reversed(range(aln2.get_alignment_length())):
    if not "-" in aln[:, pos]:
        p3p = pos
        break



aln2=aln2[:, p5p:p3p]


        
###EXTRACT SNPS and compare groups      



DICT={}
for i in range(0,len(aln2[0].seq)):
    POS=list(aln2[:, i])
    if len(set(POS))>1:
        our_clade=[]
        other_taxa=[]
        #print(POS)
        for y in range(0,len(aln2)):
            if "our_genome".upper() in str(aln2[y].id).upper():
                our_clade.append(aln2[y,i])
            else:
                other_taxa.append(aln2[y,i])
            #DICT[i][str(aln2[y].id)]=aln2[y,i]
        OUR=max({x:our_clade.count(x) for x in list(set(our_clade))}.items(), key=operator.itemgetter(1))[0]
        OTHER=max({x:other_taxa.count(x) for x in list(set(other_taxa))}.items(), key=operator.itemgetter(1))[0]
        
        if OUR!=OTHER:
            print(i,OUR,OTHER,our_clade,other_taxa)
