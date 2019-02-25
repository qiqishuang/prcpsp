import gurobipy as grb
model = grb.Model("Linear Programming")
import re
import itertools
import collections


## data ##
f = open("j301_1.txt", "r")
raw_lines  = (f.read().splitlines())
data = {}
data['A'], data['p'], data['R'], data['r'] = [], [], [], []      
line_counter = 0

for line in raw_lines:
    integer_list = list(map(int, re.findall(r'\d+', line))) # finds all the integers in a string
    if len(integer_list) == 0:
        continue
    else:
        if line_counter <= 17:
            if 'jobs' in line:
                data['n'] = integer_list[0]    
            elif 'horizon' in line:
                data['T'] = integer_list[0]    
            elif '- renewable' in line:
                data['K'] = integer_list[0]    
            line_counter = line_counter + 1
        elif int(data['n']) + 17 >= line_counter >= 18:  # PRECEDENCE RELATIONS - 49 > ... > 18
            if len(integer_list[3:]) >= 1:               # integer_list[3:] = successors
                jobnr = integer_list[0]                  # jobnr.
                data['A'].append([(jobnr, successor) for successor in integer_list[3:]])
        elif 2*int(data['n']) + 18 >= line_counter > int(data['n']) + 18:                     # 82 > ... > 50
            data['p'].append(integer_list[2])  
            data['r'].append(integer_list[3:])
        line_counter = line_counter + 1
    last_line = line_counter
    if line_counter == last_line:
        data['R'] = integer_list  
print(data)
f.close()


## Parameters ##  
def get_constants(model):
    n = int(data['n'])          ## n = 32,  dummy: 1 - 32
    T = int(data['T'])          ## T = 158
    K = data['K']               ## K = 4
    p = data['p']               ## len(p) = 32
    R = data['R']               ## len(R) = 4
    r = data['r']               ## len(r) = 32
    E = [e for e in range(n-2)] ## E = [0,...,29]
    V = [v for v in range(1,n+1)] ## V = [1,...,32]
    A = [item for sublist in data['A'] for item in sublist]
    return n, T, K, p, R, r, E, V, A
n, T, K, p, R, r, E, V, A = get_constants(model)

"""""""""""""""
initialisation
"""""""""""""""
## ok
def initial_B(A,T):
    B = [[0]*(n+1) for i in range(n+1)] ## 33 x 33 
    for i in V:
        for j in V:
            if i == j:
                B[i][j] = 0
            elif (i,j) in A: 
                B[i][j] = p[i-1]
            else:
                B[i][j] = -T
    return B

## ok
def get_F(n,A):
    if n == 2:
        F = [[i,j] for i in V for j in V for k in range(K) if i < j and (i,j) not in A and r[i-1][k] + r[j-1][k] > R[k]]
    elif n == 3:
        F_3 = [[i,j,h] for i in V for j in V for h in V if i < j < h and ((i,j) or (j,h) or (i,h)) not in A]
        F = [[i,j,h] for [i,j,h] in F_3 for k in range(K) if r[i-1][k] + r[j-1][k] + r[h-1][k] > R[k]] 
    return F

## ok
def update_B(t,B,A_0,A):               ## t: times, B: previous B, A: new A
    b = path_consistency(t-1,B,A_0,A)
    for (i,j) in A:
        b[i][j] = p [i-1]
    return b

"""""""""""""""
path consistency
"""""""""""""""
## ok
def path_consistency(t,B,A_0,A):       ## t: # times called, B: previous matrix, A: last A
    if t == 1:
        b = initial_B(A_0,T)
        count = []
        for i in V:
            for l in V:
                b[i][l] = max(B[i][j] + B[j][l] for j in V)           
                if b[i][l] != B[i][l]:
                    count.append([i,l])
        #print(len(count))
        ## if only one update            
        if len(count) == 1:
            h = count[0]
            l = count[1]
            for w in V:
                for x in V:
                     b[w][x] = max(B[w][x], B[w][h] + B[h][l] + B[l][x])

    if t > 1:
        b = update_B(t,B,A_0,A)
        count = []
        for i in V:
            for l in V:
                b[i][l] = max(B[i][j] + B[j][l] for j in V)           
                if b[i][l] != B[i][l]:
                    count.append([i,l])
        #print(len(count))
        ## if only one update            
        if len(count) == 1:
            h = count[0]
            l = count[1]
            for w in V:
                for x in V:
                     b[w][x] = max(B[w][x], B[w][h] + B[h][l] + B[l][x])
    
    return b

"""""""""""""""""""""
immediate selection
"""""""""""""""""""""
def immediate_selection(D,b,A):
    E = [(i,j) for [i,j] in D if b[i][j] >= 1 - p[j-1] and (i,j) not in A]
    if len(E) != 0:
        #print(E)
        A = list(set(E + A))
        #print(len(A))
        return A
    else:
        print('A: no update')
        return A
    
"""""""""""""""""""""
symmetric triples
"""""""""""""""""""""
def symmetric_triples(A,b,D):
    F_3 = get_F(3,A)
    ST = [[i,j] for [i,j,k] in F_3 if b[k][i] >= 1 - p[i-1] and b[k][j] >= 1 - p[j-1] and [i,j] not in D]
    if len(ST) != 0:
        print(len(ST))
        for [i,j] in ST:
                D.append([i,j])
        return D
    else:
        print('D: no update')
        return D

"""""""""""""""
edge-finding
"""""""""""""""
def frequency(D):      
    f = dict(collections.Counter(list(itertools.chain.from_iterable(D))))
    z = {}
    for k in range(2,n+1):
        z['%s'%k] = [key for key in f if f[key] >= k]   ## key = activity, where freq. of key > k
        z['%s' %k].sort()
        if len(z['%s' %k]) < k:
            del z['%s' %k]
    return z 

C = []  
def find_cliques(potential_clique=[], remaining_nodes=[], skip_nodes=[], depth=0):
    
    if len(remaining_nodes) == 0 and len(skip_nodes) == 0:
        #print('This is a clique:', potential_clique)
        C.append(potential_clique)
        return 

    for node in remaining_nodes:
        new_potential_clique = potential_clique + [node]
        new_remaining_nodes = [n for n in remaining_nodes if n in node.neighbors]
        new_skip_list = [n for n in skip_nodes if n in node.neighbors]
        find_cliques(new_potential_clique, new_remaining_nodes, new_skip_list, depth + 1)

        remaining_nodes.remove(node)
        skip_nodes.append(node)
    return C

def get_nodes(V,D):
    class Node(object):
        def __init__(self, name):
            self.name = name
            self.neighbors = []
        def __repr__(self):
            return self.name

    namespace = globals()    
    for i in V:
        namespace['N_%d' %i] = Node('%s'%i)
    for i in V:
        namespace['N_%d' %i].neighbors = [namespace['N_%d' %item[1]] for item in D if item[0] == i] + [namespace['N_%d' %item[0]] for item in D if item[1] == i]
    return

def get_clique(D):
    cliques = []
    get_nodes(V,D)
    namespace = globals()
    z = frequency(D)
    max_length_clique = int(max(z.keys()))
    for num in range(max_length_clique,1,-1): 
        clique = find_cliques(remaining_nodes=[namespace['N_%d' %i] for i in z['%s'%num]])
        cliques.append(clique)
    return cliques
            
#### find subsets ####
def subsets(C):
    subsets = set()
    for m in range(2,len(C)+1):
        subsets.update(set(itertools.combinations(C, m)))
    return subsets
######################
"ES LS"
def edge_finding(b,C):
    b_1_j, b_i_1,ind = {},{},0
    for j in C:
        print(j)
        C.remove(j)
        print(j,C)
        C_subsets = subsets(C)
        if min(b[1][i] for i in C) + sum(p[i-1] for i in C) > max(-b[i][1]+p[i-1] for i in C if i != j):
            print('True')
            b[1][j] = max(b[1][j],max((min(b[1][i] for i in C_1) + sum(p[i-1] for i in C_1)) for C_1 in C_subsets))
            b_1_j['%s' %j] = b[1][j]
            for i in C:
                if i != j:
                    b[i][1] = max(b[i][1],b[j][1]+p[i-1])
                    b_i_1['%s' %i] = b[i][1]
                else:
                    continue
        else:
            print('False')
        if 0 <= ind <= len(C)-1:
            C.insert(ind,j)
            ind +=1
        else:
            pass
    return b_1_j,b_i_1

def get_last(clique,b):
    sub_cliques = subsets(clique)
    #print(sub_cliques)
    d = {}
    for Q in sub_cliques:
        b[1][n] = max(b[1][n],min(b[1][i] for i in Q) + sum(p[i-1] for i in Q) + min(b[i][n] for i in Q))
        d['%s'%(str(Q))] = b[1][n]
    return d

####################################################################################################################
"""""""""""""""
initialisation 
"""""""""""""""
A_0 = A               ## A: 48
B = initial_B(A_0,T)  ## A: 48
D = get_F(2,A_0)      ## D: 38 (i,j) in D: j -> i  
"""""""""""""""""
local con prop
"""""""""""""""""
## 1. path consistency -> 2. immediate selection -> 'YES': update matrix B -> go to 1 - 'NO': go to 3
b = path_consistency(1,B,A_0,A)         ## 940 updates
A = immediate_selection(D,b,A_0)        #yes# A: 50
B = update_B(2,B,A_0,A)                 ## 

## 1a. - 2a.
b = path_consistency(2,B,A_0,A)         ## 866 updates
A = immediate_selection(D,b,A)          #yes# A: 54
B = update_B(3,B,A_0,A)                 ## 

## 1b. - 2b.
b = path_consistency(3,B,A_0,A)         ## 657 updates
A = immediate_selection(D,b,A)          #yes# A: 57
B = update_B(4,B,A_0,A)                 ## 

## 1c. - 2c.
b = path_consistency(4,B,A_0,A)         ## 244 updates
A = immediate_selection(D,b,A)          #no# A: 57

## 3. symmetric triples -> 'YES': go to 2 - ' NO': go to 4
D = symmetric_triples(A,b,D)            #no# D: 38
#print('D: %s'%D)

## 4. edge-finding
cliques = list(itertools.chain.from_iterable(get_clique(D)))
cliques_max = [item for item in cliques if len(item) == max(len(item) for item in cliques)]
clique_max = []
for item in cliques_max:
    item = [int(str(item[i])) for i in range(len(item))]
    item.sort()
    if item not in clique_max:
        clique_max.append(item)
print('Clique:',clique_max)


""" all False """
#d_j, d_i ={},{}
#b = path_consistency(4,B,A_0,A)  
#for C in clique_max:
#    b_1_j,b_i_1 = edge_finding(b,C)
#    if b_1_j != {}:
#        d_j['%s'%C].append(b_1_j)
#    if b_i_1 != {}:
#        d_i['%s'%C].append(b_i_1)

    
## 5. b[1][n]
last = []
for C in clique_max:
    b = path_consistency(4,B,A_0,A)         
    last.append(get_last(C,b))
#print('last: %s'%last)

m = []
for item in last:
    m.append(min(item.values()))
    b[1][n] = min(item.values())
b[1][n] = min(m)
print('b[1][n]: %s'%min(m)) ## C_max lb: 38

## earliest starting time
b_1_j, b_i_1 = {},{}
for j in V:
    b_1_j['%s' %j] = b[1][j]
print('ES: %s' %b_1_j)

## latest starting time
for i in V:
    b_i_1['%s' %i] = -b[i][1]
print('LS: %s'%b_i_1)