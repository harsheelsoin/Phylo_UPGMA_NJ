import sys
import getopt
import numpy as np
import collections


It = collections.Iterable

#Neater Output
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

#Use 'letters' object as below if verbose tree wanted printed
# letters= ['NC_012323.1_gene_3','AM489718.1_gene_16','NC_015097.1_gene_16','AP008988.1_gene_3',
# 	 'AP002929.1_gene_3','NC_002616.1_gene_3','FR751399.1_gene_16','AP008991.1_gene_3',
# 	 'AP008989.1_gene_3','AM919429.1_gene_16','AP004411.1_gene_4','AM489716.1_gene_16',
# 	 'FR751401.1_gene_16','NC_007396.1_gene_3','AB182305.1_gene_3','AP004412.1_gene_3',
# 	 'NC_010121.1_gene_4','NC_015120.1_gene_16','AP008990.1_gene_3']

#Legend:
# 'A':'NC_012323.1_gene_3','B':'AM489718.1_gene_16','C':'NC_015097.1_gene_16','D':'AP008988.1_gene_3',
# 'E':'AP002929.1_gene_3','F':'NC_002616.1_gene_3','G':'FR751399.1_gene_16','H':'AP008991.1_gene_3',
# 'I':'AP008989.1_gene_3','J':'AM919429.1_gene_16','K':'AP004411.1_gene_4','L':'AM489716.1_gene_16',
# 'M':'FR751401.1_gene_16','N':'NC_007396.1_gene_3','O':'AB182305.1_gene_3','P':'AP004412.1_gene_3',
# 'Q':'NC_010121.1_gene_4','R':'NC_015120.1_gene_16','S':'AP008990.1_gene_3'

digits = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
#digits = '0123456789'
node_dict = dict()

#Load input data from provided filename
def load_data(fn,split_lines=False):
	FH = open(fn,'r')
	data = FH.read().strip()
	FH.close()
	#print data
	if split_lines:  return data.split('\n')
	return data

def groups(L,N=3):
	R = range(0,len(L),N)
	return [L[i:i+N] for i in R]


def one_round(A,otus,count):
	print "OTUs:"
	print otus
	print "Distance matrix:"
	print A
	print
	print "Column sum:"
	div = np.sum(A,axis=1) #"divergence" matrix
	print div
	n = A.shape[0] #get number of rows
	
	# two nodes only:  we're done
	if n == 2:
		dist = A[1][0]
		nD = node_dict[otus[0]]
		nD['up'] = otus[1]
		nD['d_up'] = dist
		return None,otus

	#Find the i,j to work on using divergence:
	#based on a distance matrix relating the r taxa, calculate Q-matrix and 
	#find the pair of taxa for which Q(i,j) has its lowest value. 
	i,j = 0,0
	low_value = A[i][j]
	for r,row in enumerate(A):
		#print row 
		if r == 0:  continue
		for c,col in enumerate(row):
		#print col 
			if c >= r:  continue
			dist = A[r][c] 
			#print dist
			first = div[c]
			#print first
			second = div[r]
			#print second
			correction = (first + second)/(n-2)
			value = dist - correction
			print r, c, dist, first, second, correction, value
			if value < low_value:
				i,j,low_value = r,c,value
				
	#Merge i and j entries
	#Calculate distance of new node from tips
	new_name = digits[count]
	print
	print 'merge:', i, otus[i], j, otus[j], 'to', new_name
	
	#Dist from node[i]
	dist =  A[i][j]
	diff = div[i] - div[j]
	print 'orig dist', dist, 'div diff', diff
	dist_i = dist/2.0 + diff/(2*(n-2))
	dist_j = dist - dist_i
	#print dist_i, dist_j
	node = { 'L':otus[i], 'dL':dist_i, #"L"&"R" state for "Left" and "Ridht"
			 'R':otus[j], 'dR':dist_j } #"dL" & "dR" are left/right distances from divergence point to the left/right node. 
	node_dict[new_name] = node
	print node
	
	#Calculate distances to new node,
	# i,j assigned above
	tL = list()
	ij_dist = A[i][j]
	for k in range(len(A[0])):
		if k == i or k == j:  continue
		print 'node', otus[k], A[i][k], A[j][k], ij_dist,
		dist = (A[i][k] + A[j][k] - ij_dist)/2
		print dist
		tL.append(dist)
	print 'to new node:', tL
	
	print
	print A
	print

	#Remove columns and rows involving i or j
	if i < j:  i,j = j,i
	assert j < i
	print 'i', i, 'j', j
	sel = range(n)
	for k in [j,i]:	#Larger first
		sel.remove(k)
		print 'sel', sel
		A1 = A[sel,:]
		A2 = A1[:,sel]
		print A2
	A = A2
	print
	#Correct the otu names:
	otus = [new_name] + otus[:j] + otus[j+1:i] + otus[i+1:]
	
	#Add col at left and row at top for new node
	new_col = np.array(tL)
	new_col.shape = (n-2,1)
	A = np.hstack([new_col,A])

	new_row = np.array([0] + tL)
	new_row.shape = (1,n-1)
	A = np.vstack([new_row,A])
	print A
	print
	return A,otus
	

def build_tree(filename) :
  fn = filename 
  data = load_data(fn,split_lines=True)
  N = len(data)
  A = list()
  for line in data:
	A.append([float(n) for n in line.split()])
  otus = list(letters[:N])
  A = np.array(A)
  A.shape = (N,N)
  print A

  count = 0

  while True:
	print 'round', count
	A,otus = one_round(A,otus,count)
	if A is None:  break
	count += 1
	print A
	print
  print


  print "Printint final results:"
  print
  kL = ['L','dL','R','dR','up','d_up'] 
  for node in sorted(node_dict.keys()):
	print node, ':   ',
	for k in kL:
		nD = node_dict[node]
		if not k in nD:
			continue
		v = nD[k]
		if k in ['dL','dR']:
			v = '%3.3f' % v
		print v, ' ',
	print	

	
def main(argv) :
   #Parsing command line arguments:
   try:
	  opts, args = getopt.getopt(argv[1:],"h", ["help", "?"])
   except getopt.GetoptError:
	  print 'python nj.py <input filename> [-option [argument]]'
	  sys.exit(2)

   for opt, arg in opts:
	   if opt in ("-h", "--help", "-?"):
		   print "Reserved for help"
		   sys.exit()
   
   build_tree(args[0])
	
if __name__ == "__main__":
	main(sys.argv)