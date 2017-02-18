import collections
It = collections.Iterable

letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

def load_data(fn,do_split=False):
    FH = open(fn,'r')
    data = FH.read().strip()
    FH.close()
    if do_split:
        return data.split('\n')
    return data

def flatten(L):
    for e in L:
        if isinstance(e, It) and not isinstance(e, basestring):
            for sub in flatten(e):
                yield sub
        else:
            yield e
            
def print_list(L, n=3):
    while L:
        for e in L[:n]:
	    print e,
        print
        L = L[n:]

def list_elements(t, sort_them=True):
    return sorted(list(flatten(t)))