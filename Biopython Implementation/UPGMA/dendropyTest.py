import numpy as np 
import sys,os
import Bio
import dendropy
from dendropy.calculate import treecompare

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

aln = AlignIO.read('/home/harsheel/CPME Assignments/genealign.fasta', 'fasta')
print aln

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
print dm


constructor1 = DistanceTreeConstructor(calculator,'upgma')  #rooted=True
constructor2 = DistanceTreeConstructor(calculator,'nj')  #rooted=False

#FH=open('/home/harsheel/CPME Assignments/distanceMatrixTest.txt','r')
tree1=constructor1.build_tree(aln)
tree2=constructor2.build_tree(aln)

Bio.Phylo.draw_ascii(tree1)
#Bio.Phylo.draw(tree1)

Bio.Phylo.draw_ascii(tree2)
#Bio.Phylo.draw(tree2)

#tree = tree.as_phyloxml()
#tree1.print_newick()

print tree1
print tree2

Bio.Phylo.write(tree1, 'tree1.tre', 'newick')
Bio.Phylo.write(tree2, 'tree2.tre', 'newick')

# tree1 = Bio.Phylo.read('tree1.tre', 'newick')
# tree2 = Bio.Phylo.read('tree2.tre', 'newick')

# f1=open('tree1.tre','r')
# f2=open('tree2.tre','r')

with open('tree1.tre', 'r') as myfile:
    s1=myfile.read().replace('\n', '')

with open('tree2.tre', 'r') as myfile:
    s2=myfile.read().replace('\n', '')

print s1 
print '\n'
print s2

# establish common taxon namespace
tns = dendropy.TaxonNamespace()

# ensure all trees loaded use common namespace
t1 = dendropy.Tree.get(
        data=s1,
        schema='newick',
        taxon_namespace=tns)
t2 = dendropy.Tree.get(
        data=s2,
        schema='newick',
        taxon_namespace=tns)

## Unweighted Robinson-Foulds distance
print(treecompare.symmetric_difference(t1, t2))

# t1 = dendropy.Tree.get_from_string(tree1, 'newick')             
# t2 = dendropy.Tree.get_from_string(tree2, 'newick')

#dendropy.calculate.treecompare.robinson_foulds_distance(tree1, tree2, edge_weight_attr='length')