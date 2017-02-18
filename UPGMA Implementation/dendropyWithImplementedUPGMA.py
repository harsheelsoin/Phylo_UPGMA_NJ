import numpy as np 
import sys,os
import dendropy
from dendropy.calculate import treecompare
import Bio
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO
from upgma import run

fn = 'distanceMatrix.txt'

names2 = { 'A':'NC_012323.1_gene_3','B':'AM489718.1_gene_16','C':'NC_015097.1_gene_16','D':'AP008988.1_gene_3',
	   'E':'AP002929.1_gene_3','F':'NC_002616.1_gene_3','G':'FR751399.1_gene_16','H':'AP008991.1_gene_3',
	   'I':'AP008989.1_gene_3','J':'AM919429.1_gene_16','K':'AP004411.1_gene_4','L':'AM489716.1_gene_16',
	   'M':'FR751401.1_gene_16','N':'NC_007396.1_gene_3','O':'AB182305.1_gene_3','P':'AP004412.1_gene_3',
	   'Q':'NC_010121.1_gene_4','R':'NC_015120.1_gene_16','S':'AP008990.1_gene_3'}

node_dict = run(fn)
L = [k for k in node_dict.keys() if len(k) == 1]
otus = sorted(L)

debug = True
if debug:
    for k in node_dict:
        print k, node_dict[k]
    
for k in otus:
    tD = node_dict[k]
    if len(k) > 1:
        name = k
    if fn == 'distanceMatrix.txt':
        name = names2[k]
    else:  name = k
    tD['node_repr'] = name + ':' + '%d' % tD['up']

# build up the internal nodes by
# sorting on the length of the labels
def sort_key(s):
    return len(s)
L = sorted(node_dict.keys(), key=sort_key)
L = [k for k in L if len(k) > 1]

for k in L:
    print k
    tD = node_dict[k]
    left = node_dict[tD['left']]['node_repr']
    right = node_dict[tD['right']]['node_repr']
    try:
        up = ':%d' % tD['up']
    except KeyError:
        up = ';'
    node_repr = '(' + left + ', ' + right + ')' + up
    tD['node_repr'] = node_repr

root = ''.join(otus)
tree_text = node_dict[root]['node_repr']
#tree_text += ';'
print root, '\n'
print 'My UPGMA Implementation tree in Newick format:'
print tree_text, '\n'

aln = AlignIO.read('genealign.fasta', 'fasta')
print 'Biopython alignment read:'
print aln, '\n'

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
print 'Biopython Distance Matrix:'
print dm, '\n'

constructor1 = DistanceTreeConstructor(calculator,'nj')  #rooted=False
tree1=constructor1.build_tree(aln)

constructor2 = DistanceTreeConstructor(calculator,'upgma')  #rooted=True
tree2=constructor2.build_tree(aln)

print 'Biopython NJ Tree:'
Bio.Phylo.draw_ascii(tree1)

print '\n','Biopython UPGMA Tree:'
Bio.Phylo.draw_ascii(tree2)

#tree = tree.as_phyloxml()
#tree.print_newick()

Bio.Phylo.write(tree1, 'treeNJBiopython.tre', 'newick')
Bio.Phylo.write(tree2, 'treeUPGMABiopython.tre', 'newick')

# tree = Bio.Phylo.read('treeNJBiopython.tre', 'newick')

with open('treeNJBiopython.tre', 'r') as myfile:
    s1=myfile.read().replace('\n', '')

with open('treeUPGMABiopython.tre', 'r') as myfile:
    s2=myfile.read().replace('\n', '')

print s1, '\n'
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
t3 = dendropy.Tree.get(
        data=tree_text,
        schema='newick',
        taxon_namespace=tns)


## Unweighted Robinson-Foulds distance
print 'RF Distance between Biopython NJ tree and Biopython UPGMA Implementation tree:'
print(treecompare.symmetric_difference(t1, t2))

## Unweighted Robinson-Foulds distance
print 'RF Distance between Biopython NJ tree and my UPGMA Implementation tree:'
print(treecompare.symmetric_difference(t1, t3))

## Unweighted Robinson-Foulds distance
print 'RF Distance between Biopython UPGMA tree and my UPGMA Implementation tree:'
print(treecompare.symmetric_difference(t2, t3))

Bio.Phylo.draw(tree1)
Bio.Phylo.draw(tree2)

# t1 = dendropy.Tree.get_from_string(tree1, 'newick')             
# t2 = dendropy.Tree.get_from_string(tree2, 'newick')

#dendropy.calculate.treecompare.robinson_foulds_distance(tree1, tree2, edge_weight_attr='length')