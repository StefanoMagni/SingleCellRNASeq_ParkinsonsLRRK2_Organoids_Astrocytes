# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 17:29:58 2016

@author: tomasz.ignac
"""


import numpy as np
# import matplotlib
import matplotlib.pyplot as plt
import math
import csv

import scipy.stats as stats

from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import BernoulliNB, MultinomialNB
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA, KernelPCA
from sklearn.ensemble import RandomForestClassifier
from mpl_toolkits.mplot3d import Axes3D

#from Script_data import get_from_txt # get_gene_list,get_expr_matrix
import Script_data as sd
#import random
#import os
#import csv


def get_DR(data_matrix, g_ind=-1, d=2, off=[], t1='PCA', t2='tSNE', SaveAs1='PCA', SaveAs2='tSNE'):
    """
    Dimensionality reduction (pca, tsne, 2 or 3 dimensions)
    g_ind: index of gene used to color plots
    d: dimensions; off: indices of excluded genes (list)
    t1, t2: titles of pca and tsna plots
    """
    M = np.copy(data_matrix)
    cm = np.zeros(np.shape(M[0,:])) if g_ind==-1 else M[g_ind,:]
    M[off,:] = 0
    tsne = TSNE(n_components=d, random_state=0)   
    y_tsne = tsne.fit_transform(np.transpose(M))  
    pca = PCA(n_components=d)
    y_pca = pca.fit_transform(np.transpose(M))
    PercVarExplainPCs = pca.explained_variance_ratio_
    if d == 2:
        plt.figure()
        plt.scatter(y_pca[:,0],y_pca[:,1], c=cm, edgecolors='none')
        plt.colorbar()
        plt.title(t1)
        plt.savefig(SaveAs1)
        plt.figure()
        plt.scatter(y_tsne[:,0],y_tsne[:,1], c=cm, edgecolors='none')
        plt.title(t2)
        plt.colorbar()
        plt.savefig(SaveAs2)
    elif d == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', edgecolors='none')
        ax.scatter(y_pca[:,0], y_pca[:,1], y_pca[:,2], c=cm)
        #plt.colorbar(mappable=cm)
        plt.title(t1)
        plt.savefig(SaveAs1)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', edgecolors='none')
        ax.scatter(y_tsne[:,0], y_tsne[:,1], y_tsne[:,2], c=cm)
        #ax.colorbar()
        plt.title(t2)
        plt.savefig(SaveAs2)
    return y_pca, y_tsne, PercVarExplainPCs


def get_DR_bis(X, Y, Z=None, T=None, com=2, t1='PCA', t2='tSNE', SaveAs1='PCA', SaveAs2='tSNE'):
    M = np.hstack((X, Y, Z)) if Z is not None else np.hstack((X, Y))
    M = np.hstack((M, T)) if T is not None else M
    y_color = np.shape(X)[1]*['red'] + np.shape(Y)[1]*['blue']
    y_color = y_color + np.shape(Z)[1]*['orange'] if Z is not None else y_color
    y_color = y_color + np.shape(T)[1]*['cyan'] if T is not None else y_color
    # THE FOLLOWING TWO LINES WERE MOD BY STE AS CONFUSING AND LIKELY WRONG IN TOMAS SCRIPTS
#    y_color = np.shape(Z)[1]*['orange'] + y_color if Z is not None else y_color
#    y_color = np.shape(T)[1]*['cyan'] + y_color if T is not None else y_color
    #y_color = aux_fcol(X,Y,Z)
    tsne = TSNE(n_components=com, random_state=0)
    y_tsne = tsne.fit_transform(np.transpose(M))  
    pca = PCA(n_components=com)
    y_pca = pca.fit_transform(np.transpose(M))
    if com == 2:
        plt.figure()
        plt.scatter(y_pca[:,0],y_pca[:,1], c=y_color, edgecolors='none')
        plt.title(t1)
        plt.savefig(SaveAs1)
        plt.figure()
        plt.scatter(y_tsne[:,0],y_tsne[:,1], c=y_color, edgecolors='none')
        plt.title(t2)
        plt.savefig(SaveAs2)
    elif com == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.title(t1)
        plt.savefig(SaveAs1)
        ax.scatter(y_pca[:,0], y_pca[:,1], y_pca[:,2], c=y_color, edgecolors='none')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(y_tsne[:,0], y_tsne[:,1], y_tsne[:,2], c=y_color, edgecolors='none')
        plt.title(t2)
        plt.savefig(SaveAs2)
    return y_pca, y_tsne


#def get_DR_bis(X, Y, Z=None, T=None, com=2, t1='PCA', t2='tSNE', SaveAs1='PCA', SaveAs2='tSNE', SamplesAre4=False):
#    M = np.hstack((X, Y, Z)) if Z is not None else np.hstack((X, Y))
#    M = np.hstack((M, T)) if T is not None else M
#    if SamplesAre4 == False:
#        y_color = np.shape(X)[1]*['red'] + np.shape(Y)[1]*['blue']
#    if SamplesAre4 == True:
#        y_color = int(np.shape(X)[1]/2)*['orange'] + int(np.shape(X)[1]/2)*['cyan'] + int(np.shape(Y)[1]/2)*['red'] + int(np.shape(Y)[1]/2)*['blue']
#    y_color = y_color + np.shape(Z)[1]*['green'] if Z is not None else y_color
#    y_color = y_color + np.shape(T)[1]*['yellow'] if T is not None else y_color
#    #y_color = aux_fcol(X,Y,Z)
#    tsne = TSNE(n_components=com, random_state=0)   
#    y_tsne = tsne.fit_transform(np.transpose(M))  
#    pca = PCA(n_components=com)
#    y_pca = pca.fit_transform(np.transpose(M))
#    if com == 2:
#        plt.figure()
#        plt.scatter(y_pca[:,0],y_pca[:,1], c=y_color)
#        plt.title(t1)
#        plt.savefig(SaveAs1)
#        plt.figure()
#        plt.scatter(y_tsne[:,0],y_tsne[:,1], c=y_color)
#        plt.title(t2)
#        plt.savefig(SaveAs2)
#    elif com == 3:
#        fig = plt.figure()
#        ax = fig.add_subplot(111, projection='3d')
#        ax.scatter(y_pca[:,0], y_pca[:,1], y_pca[:,2], c=y_color)
#        fig = plt.figure()
#        ax = fig.add_subplot(111, projection='3d')
#        ax.scatter(y_tsne[:,0], y_tsne[:,1], y_tsne[:,2], c=y_color)
#    return y_pca, y_tsne


#def get_DR_tris(X, Y, Z=None, T=None, com=2, t1='PCA', t2='tSNE', SaveAs1='PCA', SaveAs2='tSNE', Color1='red', Color2='blue'):
#    M = np.hstack((X, Y, Z)) if Z is not None else np.hstack((X, Y))
#    M = np.hstack((M, T)) if T is not None else M
#    y_color = np.shape(X)[1]*[Color1] + np.shape(Y)[1]*[Color2]
#    y_color = y_color + np.shape(Z)[1]*['green'] if Z is not None else y_color
#    y_color = y_color + np.shape(T)[1]*['yellow'] if T is not None else y_color
#    #y_color = aux_fcol(X,Y,Z)
#    tsne = TSNE(n_components=com, random_state=0)   
#    y_tsne = tsne.fit_transform(np.transpose(M))  
#    pca = PCA(n_components=com)
#    y_pca = pca.fit_transform(np.transpose(M))
#    if com == 2:
#        plt.figure()
#        plt.scatter(y_pca[:,0],y_pca[:,1], c=y_color, label=('mutant','control'))
#        plt.title(t1)
#        plt.legend()
#        plt.savefig(SaveAs1)
#        plt.figure()
#        plt.scatter(y_tsne[:,0],y_tsne[:,1], c=y_color)
#        plt.title(t2)
#        plt.legend()
#        plt.savefig(SaveAs2)
#    elif com == 3:
#        fig = plt.figure()
#        ax = fig.add_subplot(111, projection='3d')
#        ax.scatter(y_pca[:,0], y_pca[:,1], y_pca[:,2], c=y_color)
#        fig = plt.figure()
#        ax = fig.add_subplot(111, projection='3d')
#        ax.scatter(y_tsne[:,0], y_tsne[:,1], y_tsne[:,2], c=y_color)
#    return y_pca, y_tsne, pca 


def get_DR_tris(X, Y, Z=None, T=None, com=2, t1='PCA', t2='tSNE', SaveAs1='PCA', SaveAs2='tSNE', Color1='red', Color2='blue', Label1=None, Label2=None):
    M = np.hstack((X, Y, Z)) if Z is not None else np.hstack((X, Y))
    M = np.hstack((M, T)) if T is not None else M
    y_color = np.shape(X)[1]*[Color1] + np.shape(Y)[1]*[Color2]
    y_color = y_color + np.shape(Z)[1]*['green'] if Z is not None else y_color
    y_color = y_color + np.shape(T)[1]*['yellow'] if T is not None else y_color
    #y_color = aux_fcol(X,Y,Z)
    tsne = TSNE(n_components=com, random_state=0)   
    y_tsne = tsne.fit_transform(np.transpose(M))  
    pca = PCA(n_components=com)
    y_pca = pca.fit_transform(np.transpose(M))
    if com == 2:
        plt.figure()
        N1 = np.shape(X)[1]
        N2 = np.shape(Y)[1]
        plt.scatter(y_pca[0:N1,0],y_pca[0:N1,1], c=y_color[0:N1], label=Label1, edgecolors='none')
        plt.scatter(y_pca[N1:N1+N2,0],y_pca[N1:N1+N2,1], c=y_color[N1:N1+N2], label=Label2, edgecolors='none')
        ax = plt.gca()
        # ax.set_facecolor('white')
        plt.title(t1)
        legend = plt.legend(frameon = 1)
        frame = legend.get_frame()
        frame.set_color('white')        
        plt.savefig(SaveAs1)
        plt.figure()
        plt.scatter(y_tsne[0:N1,0],y_tsne[0:N1,1], c=y_color[0:N1], label=Label1, edgecolors='none')
        plt.scatter(y_tsne[N1:N1+N2,0],y_tsne[N1:N1+N2,1], c=y_color[N1:N1+N2], label=Label2, edgecolors='none')
        ax = plt.gca()
        # ax.set_facecolor('white')
        legend = plt.legend(frameon = 1)
        frame = legend.get_frame()
        frame.set_color('white')  
        plt.title(t2)
        plt.savefig(SaveAs2)
    elif com == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(y_pca[:,0], y_pca[:,1], y_pca[:,2], c=y_color, edgecolors='none')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(y_tsne[:,0], y_tsne[:,1], y_tsne[:,2], c=y_color, edgecolors='none')
    return y_pca, y_tsne, pca 


def get_DR_Colored_By_User_Provided_Vector(data_matrix, ColorVecotr, MyColorMap, d=2, off=[], t1='PCA', t2='tSNE', SaveAs1='PCA', SaveAs2='tSNE', InsertColorBar=True):
    """
    Dimensionality reduction (pca, tsne, 2 or 3 dimensions)
    g_ind: index of gene used to color plots
    d: dimensions; off: indices of excluded genes (list)
    t1, t2: titles of pca and tsna plots
    """
    M = np.copy(data_matrix)
    cm = ColorVecotr
    M[off,:] = 0
    tsne = TSNE(n_components=d, random_state=0)   
    y_tsne = tsne.fit_transform(np.transpose(M))  
    pca = PCA(n_components=d)
    y_pca = pca.fit_transform(np.transpose(M))
    if d == 2:
        plt.figure()
        plt.scatter(y_pca[:,0],y_pca[:,1], c=cm, cmap=MyColorMap, edgecolors='none')
        if InsertColorBar == True:
            plt.colorbar()
        plt.title(t1)
        plt.savefig(SaveAs1)
        plt.figure()
        plt.scatter(y_tsne[:,0],y_tsne[:,1], c=cm, cmap=MyColorMap, edgecolors='none')
        plt.title(t2)
        if InsertColorBar == True:
            plt.colorbar()
        plt.savefig(SaveAs2)
    elif d == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(y_pca[:,0], y_pca[:,1], y_pca[:,2], c=cm, cmap=MyColorMap, edgecolors='none')
        plt.title(t1)
        plt.savefig(SaveAs1)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(y_tsne[:,0], y_tsne[:,1], y_tsne[:,2], c=cm, cmap=MyColorMap, edgecolors='none')     
        plt.title(t2)
        plt.savefig(SaveAs2)
    return y_pca, y_tsne


def get_probability(r1, r2=None):
    """
    This function calculate probability distribution of a vector of symbols.
    Input: list or string; or two lists r strings.(the need to be of the same length)
    Output: In case of a single input, we get probability distribution. 
    In case of two inputs, we get join distribution and marginal distribution.
    A distribution is a dictionary object. To see how it works try the following:
        P = get_probability('aaabbb')
        P
        P['a']
    """
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
    if r2 != None:
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
    if r2 == None:
        return p1
    else:
        return p1, p2, p12

    
class Compare:
    """
    This class object takes two cell populations as an input, and allows comparing them.
    A cell population is defined by a tuple: expression matrix (rows represents genes, columns are cells) and a list of genes.
    It is assumed that the order of the genes in the list is the same as in the matrix. If we use the tools for data manipulation
    this will be taken care of. 
    In order to call the class use the following: Compare((Pop1_expre_matrix, Pop1_gene_list),(Pop2_expre_matrix, Pop2_gene_list)) 
    The expression matrices can be either in the list form or numpy, gene list is a list of genes (list of strings)
    """
    def __init__(self,A,B):
        """
        Here, X, Y are expr. matrices for the two populations. x_names, y_names are corresponding gene lists
        C and c_names will be merged expr matrix and list of common genes. C contains expresions of genes that are in both X and Y
        """
        self.X = A[0] if type(A[0])==list else A[0].tolist()
        self.Y = B[0] if type(B[0])==list else B[0].tolist()
        self.C = None
        self.x_names = A[1]
        self.y_names = B[1]
        self.c_names = None
            
    def comm_genes(self):
        """produces list of genes common for the two lists"""
        aux = []
        for i in self.x_names:
            if i in self.y_names:
                aux.append(i)
        self.c_names = aux
        return self.c_names
          
    def merge(self):
        """merges the two expr. matrices into one that contains only common genes
        Additionally, we do not run com_genes before. If c_names is empty it will be run automaticaly"""
        aux=self.c_names if self.c_names!=None else self.comm_genes()
        M = []
        for i in aux:
            a = self.X[self.x_names.index(i)]
            b = self.Y[self.y_names.index(i)]
            M.append(a+b)
        self.C = M
        return np.array(M), self.c_names
    
    def merge_full(self):
        """
        Output: expression matrix with all the genes. Also includes genes expressed only in one population.
        I don't use it, we need to discuss if we need to keep it here.
        """
        M = self.C[:]
        aux = self.c_names[:]
        if M == None or aux == None:        
             M, aux = self.merge()
        x_only = []
        for i in self.x_names:
            if i not in aux:
                aux.append(i)
                x_only.append(i)
        y_only = []
        for i in self.y_names:
            if i not in aux:
                aux.append(i)
                y_only.append(i)
        x_cells = len(self.X[1])*[0]
        y_cells = len(self.Y[1])*[0]
        for i in x_only:
            a = self.X[self.x_names.index(i)]
            b = y_cells
            M.append(a+b)
        for i in y_only:
            a = x_cells
            b = self.Y[self.y_names.index(i)]
            M.append(a+b)
        return np.array(M)
         
    def MI(self, gene, num_cells_x=None, num_cells_y=None):
        """
        Mutual information (MI) between gene and variable representing populations
        gene: gene name
        num_cells_x, _y: number of cells taken from each population. By default, we take all.
        If we change this to n, top n cells will be taken.
        OUTPUT: MI score, and p-value from G test, see wikipedia
        """
        num_cells_x = len(self.X[self.x_names.index(gene)]) if num_cells_x==None else num_cells_x
        num_cells_y = len(self.Y[self.y_names.index(gene)]) if num_cells_y==None else num_cells_y
        d = num_cells_x*[0]+num_cells_y*[1]
        a = self.X[self.x_names.index(gene)][0:num_cells_x]
        b = self.Y[self.y_names.index(gene)][0:num_cells_y]
        g = a + b
        aux = []
        for i in g:
            a = 1 if i>0 else 0
            aux.append(a)
        g = aux
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
        
    def LR(self, gene, num_cells_x, num_cells_y=None):
        """as
        gjhgjhg"""
        num_cells_y = num_cells_x if num_cells_y==None else num_cells_y
        d = num_cells_x*[0]+num_cells_y*[1]
        a = self.X[self.x_names.index(gene)][0:num_cells_x]
        b = self.Y[self.y_names.index(gene)][0:num_cells_y]
        g = a + b
        x = np.transpose(np.matrix(g))
        y = np.array(d)
        B = LogisticRegression()
        B.fit(x,y)
        return B
    
    def LR_loop(self, gene_list, train, n):
        """
        Logistic regression (LR), we build a LR model based on the genes from gene_list
        train: size of the training set (value between 0 and 1). For example if we set train=0.8
        then, 80% of cells (randomly chosen) is used for training and the rest for validation of the model.
        n: number of experiments described above. Each time we chose different cells.
        Output: mean ratio of cells corretly classified 
        """
        d = len(self.X[0])*[0]+len(self.Y[0])*[1]
        g = []
        for i in gene_list:         
            a = self.X[self.x_names.index(i)]
            b = self.Y[self.y_names.index(i)]
            g.append(a+b)
        x = np.transpose(np.matrix(g))
        y = np.array(d)
        s = 0
        for i in range(n):        
            x_train, x_test, y_train, y_test = train_test_split(x, y, train_size=train)
            B = LogisticRegression()
            B.fit(x_train,y_train)
            aux = B.predict(x_test)
            aux = sum(abs(aux - y_test)) # num. of mistakes
            s = s + (np.shape(y_test)[0] - aux)/np.shape(y_test)[0]
        s = s/n
        return s
     
    def LR_any_loop(self, x, train, n, d=250*[0]+250*[1]):
        """Let's ignore it. Perhaps I'll remove it"""
        y = np.array(d)
        s = 0
        for i in range(n):        
            x_train, x_test, y_train, y_test = train_test_split(x, y, train_size=train)
            B = LogisticRegression()
            B.fit(x_train,y_train)
            aux = B.predict(x_test)
            aux = sum(abs(aux - y_test)) # num. of mistakes
            s = s + (np.shape(y_test)[0] - aux)/np.shape(y_test)[0]
        s = s/n
        return s   
    
    def NB_loop(self, gene_list, train, n, binary=True):
        d = len(self.X[0])*[0]+len(self.Y[0])*[1]
        g = []        
        for i in gene_list:        
            a = self.X[self.x_names.index(i)]
            b = self.Y[self.y_names.index(i)]
            c = a + b
            aux = []
            for j in c:
                if binary is True:
                    a = 1 if j>0 else 0
                else:
                    a = int(j) 
                aux.append(a)
            g.append(aux)
        x = np.transpose(np.matrix(g))
        y = np.array(d)
        s = 0
        for i in range(n):        
            x_train, x_test, y_train, y_test = train_test_split(x, y, train_size=train)
            B = BernoulliNB() if binary is True else MultinomialNB()
            B.fit(x_train,y_train)
            aux = B.predict(x_test)
            aux = sum(abs(aux - y_test)) #num. of mistakes
            s = s + (np.shape(y_test)[0] - aux)/np.shape(y_test)[0]
        s = s/n
        return s
        
    def RF_loop(self, gene_list, train, n, n_est=25, binary=True):
        d = len(self.X[0])*[0]+len(self.Y[0])*[1]
        g = []        
        for i in gene_list:        
            a = self.X[self.x_names.index(i)]
            b = self.Y[self.y_names.index(i)]
            c = a + b
            aux = []
            for j in c:
                if binary is True:
                    a = 1 if j>0 else 0
                else:
                    a = int(j) 
                aux.append(a)
            g.append(aux)
        x = np.transpose(np.matrix(g))
        y = np.array(d)
        s = 0
        for i in range(n):        
            x_train, x_test, y_train, y_test = train_test_split(x, y, train_size=train)
            R = RandomForestClassifier(n_estimators=n_est)
            R.fit(x_train,y_train)
            aux = R.predict(x_test)
            aux = sum(abs(aux - y_test)) #num. of mistakes
            s = s + (np.shape(y_test)[0] - aux)/np.shape(y_test)[0]
        s = s/n
        return s
    
    def plot(self, gene):
        a = self.X[self.x_names.index(gene)] + self.Y[self.y_names.index(gene)]
        plt.figure()
        plt.plot(a, 'ro')
        plt.title(gene)
        plt.xlabel('Cells ordered by population and # of transcripts')
        plt.ylabel('Transcripts')
        return 
    
    def stats(self, gene, pl=None, ret=None):
        """
        This method provides statistisc for a selected gene
        """
        a = self.X[self.x_names.index(gene)] 
        b = self.Y[self.y_names.index(gene)]
        #c = self.avg()
        mi, gt = self.MI(gene)
        mean_x = np.mean(a)
        mean_y = np.mean(b)
        an = stats.f_oneway(a,b)[1]
        kr = stats.kruskal(a,b)[1]
        lo = min(an, kr, gt)
        ratio = max(mean_y/mean_x, mean_x/mean_y) if mean_x*mean_y>0 else max(mean_x, mean_y) 
        print('P-value from anova: ', an)
        print('P-value from kr: ', kr)
        print('MI: ', mi,'with p-value= ', gt)
        print('Minimum of p-values: ',lo)
        print('mean expr. in first population: ', mean_x,'in the second:',mean_y,'ratio:',ratio)
        if pl is not None:
            self.plot(gene)
        if ret is not None:
            return an, kr, gt, lo, mi, mean_x, mean_y, ratio
        
    def avg(self):
        """
        This method calculates stats for all genes
        """
        aux=self.c_names if self.c_names!=None else self.comm_genes()

        an = [] # p-value for anova test
        kr = [] # p-value for kruskal test
        gt = [] # p-value for mutual information
        lo = [] # lowest of 3 p-values (anova, kruskal and mi)
        mi = [] # mutual information value
        mean_x = [] # mean of first population 
        mean_y = [] # mean of second population
        ratio = [] # max between mean_x/mean_y and mean_y/mean_x

        for i in aux:
            a = self.X[self.x_names.index(i)]
            b = self.Y[self.y_names.index(i)]
            aux1, aux2 = self.MI(i,len(a),len(b))
            mi.append(aux1)
            gt.append(aux2)
            mean_x.append(np.mean(a))
            mean_y.append(np.mean(b))
            if sum(a-np.mean(a))==0 or sum(b-np.mean(b))==0:
                an.append(1)
                kr.append(1)
                lo.append(1)
                ratio.append(1)
            else:
                aux3 = stats.f_oneway(a,b)[1]
                aux4 = stats.kruskal(a,b)[1]
                an.append(aux3)
                kr.append(aux4)
                lo.append(min(aux3, aux4, aux2))
                ratio.append(max(np.mean(a)/np.mean(b),np.mean(b)/np.mean(a)))
        return an, kr, gt, lo, mi, mean_x, mean_y, ratio
        
    def avg_list(self, th, stat, data=None):
        """
        data: numpy array, rows are genes
        This method creates a list of genes that have a p-value below a threshold, 
        Input: th - threshold
        stats: vector of p-values, one p-val for each gene on self.c_names
        data: additional matrix, usually with stats
        """
        if self.c_names is None or self.C is None:
            self.merge()
        lst = []
        for i in range(len(stat)):
            if stat[i] < th:
                lst.append(self.c_names[i])
        M = []
        for i in lst:
            if data is None:
                M.append(self.C[self.c_names.index(i)])
            else:
                aux = list(data[self.c_names.index(i), :])
                aux.insert(0, i)
                M.append(aux)
        return lst, M
        
   
    

        

        
"""        
St_genes, St_expr, St_expr_np = sd.get_from_txt('./DataKamil/CTR_1.txt')
Ct_genes, Ct_expr, Ct_expr_np = sd.get_from_txt('./DataKamil/ATP_1.txt')
Ct = sd.get_top_cells_matrix(Ct_expr_np,250)[0]
St = sd.get_top_cells_matrix(St_expr_np,250)[0]
Ct_norm = sd.get_normalization(Ct)[0]
St_norm = sd.get_normalization(St)[0]
A = Compare((St_norm,St_genes),(Ct_norm,Ct_genes)) 

a=A.avg()[0]
b,c=A.avg_list(a,0.01/21000)
tsne = TSNE(n_components=3, random_state=0)   
y_tsne=tsne.fit_transform(np.transpose(c))  
pca = PCA(n_components=5)
y_pca=pca.fit_transform(np.transpose(c))
y_color=250*['red']+250*['blue']
plt.figure()
plt.scatter(y_pca[:,0],y_pca[:,1], c=y_color)
plt.figure()
plt.scatter(y_tsne[:,0],y_tsne[:,1], c=y_color)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(y_pca[:,0], y_pca[:,1], y_pca[:,2], c=y_color)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(y_tsne[:,0], y_tsne[:,1], y_tsne[:,2], c=y_color)
"""


"""
#gene_c = Ct_genes.index('VIM')
#gene_s = St_genes.index('VIM')
#M[21335,:] = 0 #TAOK1 ST
#M[10639,:] = 0 #Ct, 
#M[11317,:] = 0 #St, 
#plt.title('Stimulated - PCA\npoints colored by expression of VIM')
#plt.title('Stimulated - tSNE\npoints colored by expression of VIM')
"""
