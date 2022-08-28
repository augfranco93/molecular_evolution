#------------------------------
# Molecular evolution analysis
# Augusto Franco
#------------------------------
#Library
from pandas import *
import numpy as np
from Bio import SeqIO
import os
import glob
import re,fileinput,os
#-------------------Data---------------------
genes_ids =["GENE_LIST"]

  
for i in genes_ids:  
  sequences = []    
  for seqs in SeqIO.parse("Sequence.fasta","fasta"):
    if i in seqs.id:
      sequences.append(seqs)
      SeqIO.write(sequences,"Genetic_diversity_metrics/Genes_sequences/Seqs_" + i + ".fasta","fasta")
  