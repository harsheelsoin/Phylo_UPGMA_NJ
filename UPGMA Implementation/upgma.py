import utils

def get_data_dict(fn,debug=False):
    D = dict()
    data = utils.load_data(fn,do_split=True)
    N = len(data)
    otus = utils.letters[:N]
    eL = enumerate(otus)
    L = [(o1,o2) for i,o1 in eL for o2 in otus[i+1:]]
    for o1,line in zip(otus,data):
        values = [float(s) for s in line.split()]
        for o2,n in zip((otus),values):
            k = (o1,o2)
            D[k] = n
    if debug:
        print 'starting dict'
        for k in L:
            print k, D[k]
    return D,L

def average_distance(t,D,debug=False):
    if debug:
        print 'average_distance',
        print t
    dL = list()
    t0 = utils.list_elements(t[0])
    t1 = utils.list_elements(t[1])
    if debug:
        print t0, t1,
    for u in t0:
        for v in t1:
            d = D[(u,v)]
            if debug:
                print u+v, d,
            dL.append(d)
    average = sum(dL)*1.0/len(dL)
    if debug:
        print 'avg =', average
    return average

# replace elements of t with new node
def collapse(L,node_in,debug=False):
    node = ''.join(utils.list_elements(node_in))
    if debug:
        print 'node', node
        print 'L'
        node, utils.print_list(L, n=4)
    rL = list()
    for item in L:
        left,right = item
        if debug:
            print 'left', left, 'right', right
        for tip in left:
            if tip in node:
                left = node
        for tip in right:
            if tip in node:
                right = node
        if not left == right:
            rL.append((left,right))
    rL = sorted(list(set(rL)))
    if debug:
        print 'rL after set:'
        utils.print_list(rL)
    return rL
#-----------------------------------------------
def save_node_data(t,d,node_dict,debug=False):
    # join items like (('A','B'),'C') into 'ABC'
    left = ''.join(utils.list_elements(t[0]))
    right =  ''.join(utils.list_elements(t[1]))
    new_node = ''.join(utils.list_elements(left+right))
    d *= 0.5
    node_dict[new_node] = { 'left':left, 'right':right, 'to_tips':d }
    if debug:
        print 'new node:'
        print new_node, node_dict[new_node]
    for i,child_node in enumerate([left,right]):
        if not child_node in node_dict:
            node_dict[child_node] =  {'parent':new_node, 'up':d }
        else:
            up_d = d - node_dict[child_node]['to_tips']
            node_dict[child_node].update({ 'parent':new_node, 'up':up_d} )
        if debug:
            print ['left node:','right_node:'][i]
            print child_node, node_dict[child_node]

def one_round(L,D,node_dict,debug=False):
    # find the pair of nodes with the smallest distance
    temp = [(D[k],k) for k in L]
    d,t = sorted(temp)[0]
    if debug:
        print 'closest:'
        print t,d
    save_node_data(t,d,node_dict,debug=debug)
    # collapse the list
    rL = sorted(collapse(L,t))
    if debug:
        print 'elements after joining:'
        utils.print_list(rL,n=1)
    # calculate cluster distances
    for t in rL:
        if not t in D:
            D[t] = average_distance(t,D,debug=debug)
    return rL,t

def run(fn,debug=False):
    D,L = get_data_dict(fn,debug=debug)
    node_dict = dict()
    counter = 0
    while L:
        counter += 1
        if debug:
            print '\nround', counter
        L,t = one_round(L,D,node_dict,debug=debug)
    if debug:
        print t, '\n'
        for k in node_dict:
            print k, node_dict[k]
    return node_dict    
#-----------------------------------------------
if __name__ == '__main__':
    
    fn = 'distanceMatrix.txt' 	
    node_dict = run(fn,debug=True)
    
    print '***'
    def f(k):  return len(k)
    L = sorted(node_dict.keys(),key=f)
    for k in L[8:]:
        print k,
    print
