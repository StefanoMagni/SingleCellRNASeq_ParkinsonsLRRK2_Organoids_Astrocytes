# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 17:29:58 2016

@author: tomasz.ignac
"""


import numpy as np
import math
#from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2
#import random

 
def get_probability(r1, r2=None):
    r1 = r1.tolist() if type(r1)==np.ndarray else r1 
    r2 = r2.tolist() if type(r2)==np.ndarray else r2 
    aux1 = set(r1)
    p1 = []
    for i in aux1:
        p1.append([i, r1.count(i)/len(r1)])
    p1 = dict(p1)    
    if r2!=None:
        p2 = []
        aux2 = set(r2)
        for i in aux2:
            p2.append([i, r2.count(i)/len(r2)])
    else:
        p2 = None 
    if r2!=None:
        aux12 = []
        for i in aux1:
            for j in aux2:
                aux12.append((i, j))
        mer = []         
        for i, iaux in enumerate(r1):         
            mer.append((r1[i],r2[i]))
        p12 = []
        for i in aux12:
            p12.append([i, mer.count(i)/len(r1)])
    else:
        p12 = None
    p2 = dict(p2) if p2!=None else None
    p12 = dict(p12) if p12!=None else None
    if r2==None:
        return p1
    else:
        return p1, p2, p12

def MI(d, g=None):
    g = d if g is None else g
    a,b,c = get_probability(d,g)
    aux1 = 0
    for i,j in c.items():
        aux1 = aux1 + j*math.log(j,2) if j>0 else aux1
    aux2 = 0
    for i,j in a.items():
        for k,l in b.items():
            aux2 = aux2 + j*l*math.log(j*l,2) if j*l>0 else aux2
    mi = aux1 - aux2
    aux = mi/math.log(math.exp(1),2)#change the base for the G-test 
    gt = 1-chi2.cdf(2*len(d)*aux,((len(a)-1)*(len(b)-1)))
    return mi, gt
        
def ID(a, b, c):
    out = 0
    return out
       
