
# coding: utf-8

# In[1]:

from __future__ import division
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


# In[14]:

def a_run_workers(workers,deadline,x,z,zdot,w,theta,penalty,x0,u,v,C, dis, alpha, Av, Pr, Rj,Jp, Bound, optimalvalue):
    
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
        
        if alpha  ==1:
            fairzdot=  np.sum(w*np.log(np.sum(zdot, axis = 1)))
            
        else:
            
            fairzdot = np.sum(w*np.sum(zdot, axis = 1)**(1-alpha)/(1-alpha))
        FF[timeproj] = fairzdot
        FM[timeproj]= fairzmin
        Gaps[timeproj] = np.abs(fairzdot - optimalvalue)/optimalvalue *100

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


# In[15]:

def re_init(*l):
    for a in l:
        a-=a


# In[ ]:

def run_one(name,deadline):
    Hst,Gps,It = {},{},{}
    WorkersList,x,z,zdot,w,theta,Penalties,x0,u,v,C, dis, alpha_list, Av, Pr, Rj,Jp, Bounds, info = load_obj('Instance_tree'+str(name))
    OptimalValue = load_obj('Opt_Tree'+str(name))
    for alpha in alpha_list:
        for dd in [Hst,Gps,It]:
            dd[alpha] = []
        optimal = OptimalValue[alpha][max(OptimalValue[alpha].keys())]
        print 'Doing tests on instance '+str(name)+'for alpha = '+str(alpha)
        for (n_cores, workers) in WorkersList:
            
            print 'Now operating for '+str(n_cores)+' cores...'
    
            N,x,zdot, zmin,gap,N,FF,FM,Zs,Zmins,Gaps =a_run_workers(workers,
                                                                    deadline,
                                                                    x,z,zdot,w,theta,
                                                                    Penalties[alpha],
                                                                    x0,u,v,C, dis,
                                                                    alpha, Av, Pr, Rj,Jp,
                                                                    Bounds[alpha], optimal)
            Hst[alpha].append((n_cores,Gaps))
            Gps[alpha].append((n_cores,gap))
            It[alpha].append((n_cores,N))
            re_init(zdot,x,z,u,v)
            print 'Finished the core '+str(n_cores)
            save_obj((Gps,Hst,It), 'Tree_Result'+str(name))
        print 'Finished the alpha '+str(alpha)

