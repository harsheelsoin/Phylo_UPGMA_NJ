import numpy as np 
import sys,os
import dendropy
from dendropy.calculate import treecompare
import Bio
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO


aln = AlignIO.read('genealign.fasta', 'fasta')
print 'Biopython Alignment Read:'
print aln, '\n'

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
print 'Biopython Distance Matrix:'
print dm, '\n'

constructor = DistanceTreeConstructor(calculator,'nj')  #rooted=False
tree=constructor.build_tree(aln)

print 'Bipython NJ Tree Plot:'
Bio.Phylo.draw_ascii(tree)
print tree

Bio.Phylo.write(tree, 'treeNJBiopython.tre', 'newick')

# tree = Bio.Phylo.read('treeNJBiopython.tre', 'newick')

with open('treeNJBiopython.tre', 'r') as myfile:
    s1=myfile.read().replace('\n', '')

#Reading my NJ Tree from text file and storing in python oject 's2'
with open('newickNJTreeFromDistMatrix.txt', 'r') as myfile:
    s2=myfile.read().replace('\n', '')

print 'Biopython NJ Tree in Newick format:'
print s1, '\n' 
print 'My NJ tree in Newick format:'
print s2, '\n'

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
print 'Distance between Biopython NJ tree and my NJ Implementation tree:'
print(treecompare.symmetric_difference(t1, t2))

Bio.Phylo.draw(tree)