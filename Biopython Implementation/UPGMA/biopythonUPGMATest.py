import numpy as np 
import sys,os
import Bio
import dendropy

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO

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
#text_file = open("distanceMatrix.txt", "w")

for seq1 in seqList:
	mutationsList=[]
	for seq2 in seqList:
		mutations=0
		for i,j in zip(seq1,seq2):
			if i != j and i != '-' and j != '-':
				mutations+=1
		mutationsList.append(mutations)
        	#text_file.write("%3d " % mutations)
	#text_file.write('\n')
    	distanceMatrix.append(mutationsList)

print distanceMatrix

aln = AlignIO.read('genealign.fasta', 'fasta')
print aln

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
print dm

#constructor = DistanceTreeConstructor(calculator,'nj')  #rooted=False
constructor = DistanceTreeConstructor(calculator,'upgma')  #rooted=True

#FH=open('/home/harsheel/CPME Assignments/distanceMatrixTest.txt','r')
tree=constructor.build_tree(aln)
Bio.Phylo.draw_ascii(tree)
Bio.Phylo.draw(tree)

print tree
