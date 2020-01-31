# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 22:10:54 2020

@author: lliu01
"""

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Alphabet import generic_rna
try:
    f=open('C:/Users/LLIU01/OneDrive - Danaher/Li/1-Projects/7- BioInformatics/dna2.fasta')
except IOError:
    print("File does't exist")

line_count = 0
for line in f:
    if line[0]=='>':
        line_count +=1
print(line_count)

file = "C:/Users/LLIU01/OneDrive - Danaher/Li/1-Projects/7- BioInformatics/dna2.fasta"
 
from Bio import SeqIO
for seq_record in SeqIO.parse(file, "fasta"):
    print(len(seq_record))
    print(seq_record.id)

seq_id_list = []
seq_len_list = []
for seq_record in SeqIO.parse(file, "fasta"):
    seq_id = seq_record.id
    seq_id_list.append(seq_id)
    seq_len = len(seq_record)
    seq_len_list.append(seq_len)
#print(seq_id_list)
#print(seq_len_list)

ma = max(seq_len_list)
max_index = [i for i, j in enumerate(seq_len_list) if j == ma]
mi = min(seq_len_list)
min_index = [i for i, j in enumerate(seq_len_list) if j == mi]
for i in range(len(max_index)):
    print(seq_id_list[max_index[i]], seq_len_list[max_index[i]])
for i in range(len(min_index)):
    print(seq_id_list[min_index[i]], seq_len_list[min_index[i]])
  
'''
seq_list.append(seq_record)
print(str(seq_list[1].seq))    
for seq in seq_list:
    seqs = str(seq)
seq1 = seq_list[0]
print(str(seq1.seq))

st = 'ATCACCTTATGCCTTTAATAGTGAGGGGC'
substr = slice(1,len(st),3)
print(str[substr])
'''

file = "C:/Users/LLIU01/OneDrive - Danaher/Li/1-Projects/7- BioInformatics/dna2.fasta"
def findORF(myfile):
    for seq_record in SeqIO.parse(myfile, "fasta"):
        seqs = str(seq_record.seq)
        print(seqs)
        start = -1
        stop = -1
        for i in range(1, len(seqs)):
            frame = st[i:i+3]
            if frame == 'ATG':
                start = i+1
                break
            i = i+3
        for i in range(start+2, len(seqs)):
            frame = st[i:i+3]
            if frame == 'TAA' or frame == 'TAG' or frame == 'TGA':
                stop = i+1
                break
            i = i+3
        print(start, stop)
findORF(file)

st = 'GCGGTCCCGGCGCCGCAGGCCGTCCGGCTCCTGCAGCGCGCCGAACCGGGTCTCGCGGTGATTGCCCAGCGTACCGAGATGCGCCCGGCCTGGGCCGTGATGGCGCAGTGCGGCC'

for i in range(1,len(st)):
    frame = st[i:i+3]
    if frame == 'ATG':
        start = i+1
        print(start)
        break
    i = i+3
for i in range(start+2,len(st)):
    frame = st[i:i+3]
    if frame == 'TAA' or frame == 'TAG' or frame == 'TGA':
        stop = i+1
        print(stop)
        break
    i = i+3
