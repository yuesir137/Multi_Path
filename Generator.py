
# coding: utf-8
from __future__ import division
# In[1]:

import csv
import networkx as nx
import matplotlib.pyplot as plt

import os
import tempfile
import shutil
import scipy.stats as stats
import math
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import time
from joblib import Parallel, delayed
import random as rd
import pickle

def save_obj(obj,name):
	with open('datas/'+name + '.pkl', 'wb') as f:
        	pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('datas/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


# In[142]:

def generate_source(G,u,v, max_split,x,y,z,Paths, out):
    L = G.nodes()
    P = []
    number = 0
    qq = [p for p in nx.all_shortest_paths(G,u,v)]
    nb_chemins = len(qq)
    p = qq.pop()
    P.append(p)
    for number in range(1,max_split):
        
        #for j in G[u]:
        #    G[u][j]['weight'] = rd.uniform(1,10)
        #print G[12][131]['weight']
            #G[j][i]['weight'] = rd.uniform(0,100)
        
        #p  = nx.shortest_path(G,u,v,'weight')
        if qq==[]:
            break
        p = qq.pop(rd.randint(0,len(qq)-1))
        
        if out==1:
            
            p.pop(-1)
            
        elif out ==2:
            
            p.pop()
            
        
        if p not in P:
            P.append(p)
    if x.keys()==[]:
        xs = 0
    else:
        xs = max(x.keys())+1
    x[xs] = 0
    Paths[xs] = 0
    for rr in range(0,len(P)):
        y[xs,rr] = 0
        p = P[rr]
        Paths[xs]+=1
        U = p.pop(0)
        while p != []:
            V = p.pop(0)
            if (U,V) in z:
                z[U,V][xs,rr] = 0
            else:
                z[U,V] = {(xs,rr):0}
            U = V
        longueur = len([(U,V) for (U,V) in z if (xs,rr) in z[U,V]])
        if longueur ==0:
            
            print xs,rr, 'uses', longueur , 'links', P[rr], 'although we had', nb_chemins,'paths', 'from',u, 'to',v
    #return x,y,z,P
        
            
        


# In[115]:

g=nx.Graph()
k = 16 #can be k = 8 (tiny) or k = 32 (huge)
fig = plt.figure(figsize = (20,12))
#extract and add AGE_GROUP nodes in graph
f1 = csv.reader(open("FAT/fat_tree-connections-K"+str(k)+".txt","rb"))
for row in f1: 
    g.add_edges_from([(row[0],row[1])])
#colors = ['pink' for kk in range(0,k**3//4)]+['green' for kk in range(k**3//4,k**3//4+k**2//2)]+['red' for kk in range(
#        k**3//4+k**2//2,                                                                                            
#        k**3//4+k**2//2+k**2//2)] +["blue" for kk in range(
#        k**3//4+k**2//2+k**2//2,int(len(g.nodes())))]
h = nx.Graph()
h.add_edges_from([(int(u),int(v)) for (u,v) in g.edges()])
ROOT = max(int(u) for u in g.nodes())+1
print 'ROOT=', ROOT

Servers = range(0,k**3//4)
Edges = set()
Aggregation = set()
Core = set()
for u in Servers:
    for v in h.neighbors(u):
        Edges.add(v)
for u in Edges:
    for v in h.neighbors(u):
        if v not in Servers:
            Aggregation.add(v)
for u in Aggregation:
    for v in h.neighbors(u):
        if v not in Edges:
            Core.add(v)
h.add_edges_from([(u,ROOT) for u in Core])
#colors = ['purple' for u in Servers]+['blue' for u in Edges]+['green' for u in Aggregation]+['orange' for u in Core]            
#colors.append('red')
#nx.draw(h,node_list = sorted(h.nodes()),node_color=colors)


# In[131]:

def generate_index(x,y,z,Paths):
    Pr = {}
    Rj = {}
    Jp = {}
    for j in z:
        Rj[j] = z[j].keys()
    for r in x:
        Pr[r] = range(0, Paths[r])
    nq = max(Paths[r] for r in Paths)
    for r in x:
        for p in Pr[r]:
            Jp[r,p] = []
            for j in z:
                if (r,p) in z[j]:
                    Jp[r,p].append(j)
    Av = np.zeros((len(x), nq))
    for (r,p) in Jp:
        Av[r][p] = len(Jp[r,p])+1
        if Jp[r,p] == []:
            print 'WARNING', r,p, 'is None'
    return Pr, Rj, Jp, Av


# In[6]:

def a_instantiate(Pr,Rj,Jp):
    nr = max(len(Pr[r]) for r in Pr)
    x = np.zeros((len(Pr), nr))
    u = np.zeros(x.shape)
    z = np.zeros((len(Rj),len(Pr),nr))
    v = np.zeros(z.shape)
    zdot = np.zeros(x.shape)
    theta = np.zeros(x.shape)
    x0 = np.zeros(x.shape)
    Dis = np.zeros(x.shape)
    return x,u,z,v,zdot,theta,x0,Dis


# In[62]:

def compute_bound(alpha,w,Jp,zdot,C, Pr,Rj):
    
    Bound = {}
    LL = {}
    minC = {}
    
    
    for (r,p) in Jp:
        if Jp[r,p] ==[]:
            minC[r,p] = 0
        else:
            minC[r,p] = min([C[jj] for jj in Jp[r,p]])
    
    a = {}
    for r in range(zdot.shape[0]):
        a[r] = sum(minC[r,p] for p in Pr[r])
    www = np.sum(w)
    if True:
        for r in Pr:
            for p in Pr[r]:
                l = set()
                for j in Jp[r,p]:
                    for (s,q) in Rj[j]:
                        l.add(s)
            l.add(r)
            LL[r] = l
            
    if alpha ==1:
        for rr in Pr:
            #print 'L[rr]=', LL[rr]
            W = sum(w[ss] for ss in LL[rr] )
            #print 'W=', W
            Bound[rr] = float(w[rr])/W * a[rr]
    elif alpha>1:
        U = (min([w[ss]*a[ss]/sum(w[tt] for tt in LL[ss]) for ss in Pr]))
        #print U
        for rr in Pr:
            #U = (min([w[ss]*a[ss]/sum(w[tt] for tt in LL[ss]) for ss in LL[rr]]))
            Bound[rr] = U**(1-1/alpha) * (w[rr]*a[rr]/sum(w[tt] for tt in LL[rr]))**(1/alpha)
    elif 0<alpha<1:
        Bmaxalpha = (max(a[rr] for rr in a))**(1-1/alpha)
        for rr in Pr:
            Bound[rr] =(w[rr]*a[rr]/sum(a[ss]**(1-alpha)*w[ss] for ss in LL[rr]))**(1/alpha)
            
    clam =1./alpha * np.sqrt(1./(min(w[r]/a[r]**(alpha+1) for r in Pr )*max(w[r]/(Bound[r])**(alpha+1) for r in Pr)))    
        
    return Bound, clam


# In[8]:

def convert(Rj,Jp):    
    jn = 0
    Rjj = {}
    for j in Rj:
        Rjj[jn] = Rj[j]
        jn+=1
    Jpp = {}
    for (r,p) in Jp:
        Jpp[r,p] = []
        for j in Rjj:
            if (r,p) in Rjj[j]:
                Jpp[r,p].append(j)
    return Rjj,Jpp


# In[79]:

def assign_to_workers(Rj,Pr,number):
    	workers = {}
   	lenpr = len(Pr)
    	j_modulo = len(Rj)%number
    	j_quo = int(len(Rj)/number)
    	for i in range(number-1):
        	workers[i] ={'proj':[j for j in range(i*j_quo,(i+1)*j_quo)], 'time' :0.0, 'prox':[]}
    	workers[number-1] = {'time':0, 'proj':[j for j in range((number-1)*j_quo,len(Rj))], 'prox':[]}
    	Total = Pr.keys()
	Already = set([])
	liste = workers.keys()
	for r in Total:
		rd.shuffle(liste)
		for d in liste:
			for j in workers[d]['proj']:
				for (s,p) in Rj[j]:
					if ((s ==r) and (r not in Already)):
						
						workers[d]['prox'].append(r)
						Already.add(r)
						break
	print 'assigned a total of', len(Already), 'requests. Hope its fine... I had', len(Pr), 'users in the beginning'
	if len(Already) != len(Pr):
		for r in Pr:
			if r not in Already:
				print r, 'is missing'
				for p in Pr[r]:
					print p, [j for j in Rj if (r,p) in Rj[j]]
	return workers, len(Already), len(Pr)


# In[144]:

def generate_problem(graph,graph_type = 8, n_requests = 1, n_split=  1, alpha_list = [1,2,4]):
    x,y,z,Paths  = {},{},{},{}
    number =0
    Bounds = {}
    n_servers = graph_type**3//4
    root = len(graph.nodes()) -1
    for u in range(n_servers):
        generate_source(graph,u,root, n_split,x,y,z,Paths,1)
        generate_source(graph,u,root, n_split,x,y,z,Paths,2)
        v = rd.randint(0,n_servers -1)
        while v==u:
            v = rd.randint(0,n_servers -1)
        if number < n_requests:
            
            generate_source(graph,u,v, n_split,x,y,z,Paths, 0)
            number +=1
    Pr, Rj, Jp, Av=generate_index(x,y,z,Paths)
    x,u,z,v,zdot,theta,x0,dis = a_instantiate(Pr,Rj,Jp)
    w = np.full(x.shape[0],1)
    C = np.full(len(Rj),100)
    Rj,Jp = convert(Rj,Jp)
    Penalties ={}
    for alpha in alpha_list:
        Bounds[alpha],clam =compute_bound(alpha,w,Jp,zdot,C, Pr,Rj)
        Penalties[alpha] = clam
    return Pr,Rj,Jp,Av,x,u,z,v,zdot,theta,x0,dis,w,C,Penalties,Bounds


# In[20]:

def solve_real_polynomial(a):    
    
    b = np.roots(a)
    c = np.copy(b)
    
    b = b[b.imag ==0].real
    return b[b>0]
    


def proximity_operator(r,alpha,w,x0,x1,u, theta, penalty, Bound):
    coeffs = np.zeros(alpha+2)
    coeffs[0] = 1
    coeffs[-1] = -len(x0)*penalty*w
    coeffs[1] = -sum(u)
    y =max(0,solve_real_polynomial(coeffs)[0])**alpha
    sol = u+ penalty*w/y
    #print 'solution',r,y,sol
    return sol


# In[12]:

def a_simplex_projection(C, Dis, zin,L):
    """Exact projection onto a simplex - L is the list of index names of vector z"""
    
    L.sort()
    z = np.copy(zin)
    r0,p0 = z.shape
    
    j = np.inf
    for i in (L):
        if j==i[0]:
            pass
        else:  
            j = i[0]
            i = i[0]
            z[i][:] = np.maximum(Dis[i], z[i])
    if (np.sum(z[r][p] for (r,p) in L) > C):
        j = np.inf
        for i in (L):
            if j==i[0]:
                pass
            else:
                
                j = i[0]
                i = i[0]
                z[i][:] = np.maximum(0, z[i] - Dis[i])     
        CC = C - np.sum(Dis[r][p] for (r,p) in L)
        u =sorted(L, key = lambda rr: z[rr[0]][rr[1]], reverse = True)
        b = False
        lenU = len(u)
        q=0
        value = np.inf
        
        while ((q<lenU) & (value >0)):
            q+=1
            value =(z[u[min(lenU-1,q)][0]][u[min(lenU-1,q)][1]] + 1./(q+1)*(CC -sum(z[u[min(lenU-1,k)][0]]
                                                                                    [u[min(lenU-1,k)][1]] 
                                                                                        
                                                                                    for k in range(0,q+1))))
        if q==0:
            q+=1
        rho = 1./(max(1,q)) * (CC - sum(z[u[k][0]][u[k][1]] for k in range(0,q)))
        #print rho
        j = np.inf
        for i in sorted(L):
            if j==i[0]:
                pass
            else:
                j = i[0]
                i = i[0]
                z[i][:] = np.maximum(z[i]+rho,0) +Dis[i]
                
        
        
    return z


# In[166]:

def create_problem(graph,graph_type = 8, n_requests = 1, n_split=  1, alpha_list = [1,2,4], n_cores = [1], name = 1):
    Pr,Rj,Jp,Av,x,u,z,v,zdot,theta,x0,dis,w,C,Penalties,Bounds= generate_problem(graph,graph_type,
                                                                                n_requests,
                                                                                n_split,
                                                                                alpha_list,
                                                                               )
    WorkersList = [] 
    info = "WorkersList,x,z,zdot,w,theta,Penalties,x0,u,v,C, dis, alpha_list, Av, Pr, Rj,Jp, Bound, info"
    for cores in n_cores:
        workers,aa,bb = assign_to_workers(Rj,Pr,cores)
        WorkersList.append((cores, workers))
    
    save_obj((WorkersList,x,z,zdot,w,theta,Penalties,x0,u,v,C, dis, alpha_list, Av, Pr, Rj,Jp, Bounds, info), 'Instance_tree'+str(name))
    return aa,bb
    


# In[14]:

def a_run_workers(workers,deadline,x,z,zdot,w,theta,penalty,x0,u,v,C, dis, alpha, Av, Pr, Rj,Jp, Bound):
    
    proxu = np.zeros(u.shape)
    Zs ={}
    Zmins = {}
    FF = {}
    FM = {}
    Gaps = {}
    zmin = np.copy(zdot)
    parallel_time = 0
    consensus_time = 0
    proxv = np.zeros(v.shape)
    nj,nr,nq = z.shape
    gap =100
    for r in Pr:
        for p in Pr[r]:
            Av[r][p] = len(Jp[r,p]) +1
    N = 0
    zdotold = np.copy(zdot)
    timeprox = 0
    timeproj = 0
    Start = time.time()
    while (timeproj <deadline):
        start = time.time()
        zdotold[:] = zdot
        proxu = zdot - penalty*u
        proxv = np.full((nj,nr,nq),zdot) - penalty*v 
        for d in workers:
            start = time.time()
            for r in workers[d]['prox']:
                x[r][0:len(Pr[r])] = proximity_operator(r,alpha,w[r],x0[r][0:len(Pr[r])],x0[r][0:len(Pr[r])],
                                                    proxu[r][0:len(Pr[r])], theta[r][0:len(Pr[r])], penalty, Bound[r])
        
            

           
            for j in workers[d]['proj']:
                z[j] = a_simplex_projection(C[j], dis,proxv[j],Rj[j])
            end = time.time()
            timeprox = max(timeprox,end-start)
        parallel_time +=timeprox
        #print 'workers worked for', timeprox,'seconds.'
        Zs[N] = np.zeros(zdot.shape)  
        start = time.time()
        for r in Pr:
            for p in Pr[r]:
                value = 1/(Av[r][p])*(sum(z[j][r][p]  for j in Jp[r,p]) +x[r][p] )
                #zdot[r][p] = 1/(Av[r][p])*(sum(proxv[j][r,p] +penalty*v[j][r][p] for j in Jp[r,p]) +x[r][p] +penalty* u[r][p])
                zdot[r][p] = value
                Zs[N][r][p] = value
                
        
        #print zdot[0][0],
        u +=1./penalty*(x-zdot)
        v +=1./penalty*(z - np.full((nj,nr,nq),zdot))
        end = time.time()
        timeproj +=end-start
        consensus_time+=end-start
        #print 'algorithm been done overall for', timeproj,'seconds'
        timeproj+=timeprox
        timeprox =0
        zmin = np.min(z,axis = 0)
        Zmins[N] = zmin
        N+=1
        zmin = np.min(z,axis = 0)
        
        gap = np.linalg.norm(np.abs(np.sum(x,axis=1) -np.sum(zdot,axis=1)),2)+np.linalg.norm(np.sum(zdot,axis= 1) - np.sum(zdotold,axis =1),2)
        Gaps[timeproj] = gap
        if alpha  ==1:
            fairzdot=  np.sum(w*np.log(np.sum(zdot, axis = 1)))
            fairzmin=  np.sum(w*np.log(np.sum(zmin, axis = 1)))
        else:
            fairzmin = np.sum(w*np.sum(zmin, axis = 1)**(1-alpha)/(1-alpha))
            fairzdot = np.sum(w*np.sum(zdot, axis = 1)**(1-alpha)/(1-alpha))
        FF[timeproj] = fairzdot
        FM[timeproj]= fairzmin
        #penalty *=rd.uniform(.9,1.1)
        
        #print penalty,
    print
    print '2-norm between consensus and zmin=',gap
    print 'Times:'
    print 'Proximal operators =', parallel_time
    print 'Consensus =', consensus_time
    print 'Estimated parallel total =', timeproj
    print 'Total =', time.time() - Start
    print 'Iterations:'
    print N
    return N,x,zdot, zmin,gap,N,FF,FM,Zs,Zmins,Gaps


# In[189]:


if  __name__ =='__main__':
    number_of_instances = 5
    for i in range(number_of_instances):
        create_problem(graph = h,graph_type = 16, n_requests = 100000, n_split=  4, alpha_list = [1,2,4], n_cores = [1,2,4,16,32], name = i)


# In[15]:

def re_init(*l):
    for a in l:
        a-=a


# In[ ]:



