from rpy2 import robjects
from rpy2.robjects.packages import importr
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
print tree_text

ape = importr('ape')
grdevices = importr('grDevices')
tr = ape.read_tree(text=tree_text)
pL = fn.split('.')
ofn = '.'.join(pL[:-1]) + '.pdf'

grdevices.pdf(ofn)
# for distanceMatrix data, manually adjusted x-axis
ape.plot_phylo(tr,cex=2,edge_width=2)
#ape.plot_phylo(tr,cex=2,edge_width=2,x_lim=100)
ape.axisPhylo()
grdevices.dev_off()