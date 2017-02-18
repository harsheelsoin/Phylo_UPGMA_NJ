import numpy as np
import os

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

fp=open("genealign.fasta",'r')
seqList=[]
for name, seq in read_fasta(fp):
        #print(name, seq)
        seqList.append(seq)

distanceMatrix=[]
text_file = open("distanceMatrix.txt", "w")

for seq1 in seqList:
	mutationsList=[]
	for seq2 in seqList:
		mutations=0
		for i,j in zip(seq1,seq2):
			if i != j and i != '-' and j != '-':
				mutations+=1
		mutationsList.append(mutations)
        	text_file.write("%3d " % mutations)
	text_file.write('\n')
    	distanceMatrix.append(mutationsList)

print distanceMatrix
#print len(distanceMatrix)
