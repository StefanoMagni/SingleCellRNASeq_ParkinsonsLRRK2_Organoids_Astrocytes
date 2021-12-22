# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 17:29:58 2016

@author: tomasz.ignac
"""

"""
Here we have functions to manipulate data. 
They are used to read the data, perform the cleaning and save output
"""



import numpy as np
import math
import csv


def get_from_csv(name, deli = ','):
    """ 
    delimiter is set to ',' but it can be changed by the user 
    the function is suposed to read only matrices of number, no text, e.g., gene names
    """
    M = []
    with open(name) as csvfile:
        readCSV = csv.reader(csvfile, delimiter = deli)
        for i in readCSV:
            r = []
            for j in i:
                r.append(float(j))
            M.append(r)
    return np.array(M)

def get_from_txt(name):
    """
    reads data from text files delivered by Drop_Seq bioinf pipline
    first column contains gene names
    Input: name=file
    Output: gene list, expression matrix in list for, expression matrix in numpy format
    """
    f = open(name)
    M = f.readlines()
    f.close()
    G = get_gene_list(M) #described below
    E, E_np = get_expr_matrix(M) #described below
    return G,E,E_np 
    
def get_gene_list(M):
    """
    Input: text file delivered by Drop_Seq bioinf pipline, where
    the first column contains gene names
    Output: list of gene names from the first column
    """
    gn = [] #gn stands for gene names
    for line in M[1:]:
        name = ''
        for i in line:
            if i==',' or i=='\t' or i=='\n':
                break
            name = name + i
        gn.append(name)
    return gn
    
def get_expr_matrix(raw):
    """
    Input: text file delivered by Drop_Seq bioinf pipline, where
    the first column contains gene names
    Action: for each row(gene) calls function that extracts the actuall expression acros cells
    Output: expr. matrix in list form, expr. matrix in numpy form
    """
    matrix = []
    for i in raw[1:]:
        matrix.append(get_gene_values(i)) #described below
    return matrix, np.array(matrix)

def get_gene_values(raw, name=None):
    """
    Input: row of a matrix where the first column is gene name
    Output: gene expressions across all cells (cols.) in list form
    """
    v = []   
    aux = ''
    for i in raw:
        if i==',' or i=='\t' or i=='\n':
            try:
                v.append(int(aux))
            except ValueError:
                pass
            aux = ''
        else:
            aux = aux + i     
    return v
       
def get_binary(M):
    """
    Simple function that performs binarization based on the median of gene expression
    Input: matrix where rows are genes, columns represent cells
    Output: corresponding matrix with 0/1 expressions
    """
    M = np.array(M)
    aux = np.median(M,axis=1)
    M_bin = np.empty([np.shape(M)[0],np.shape(M)[1]])
    for i in range(np.shape(M)[0]):
        for j in range(np.shape(M)[1]):
            M_bin[i,j] = 1 if M[i,j]>aux[i] else 0
    return np.array(M_bin)

def get_knee_plot(expr):
    """
    This function is used to estimate quality of data, i.e., diftribution of trancript across cells
    I plot outcomes of the function to decide how many cells pick for the further analysis
    Input: numpy expression matrix
    Action: the cells are ordered by the total expression (sum of all reads, i.e., columns are ordered by their sums).
    Output: 1. Cumulative expression, vector of sum of all reads for each cell 
    """
    c = []
    d = []
    for i in list(range(expr.shape[1])):
        c.append(sum(expr[:,i]))
    c = sorted(c, reverse=True)
    aux = sum(c)
    for i in list(range(len(c))):
        d.append(sum(c[:i])/aux) #here, we calculate the cumulative expr
    return  d, c
    
def get_top_cells_matrix(M,n):
    """
    Input: np.array from get_expr_matrix, where rows are genes, columns represent cells; number of cells
    Action: we select n columns(cells) with the highest number of reads (column sum)
    Output: A matrix with n columns ordered by the total number of reads. 
            There are two outputs: np.array and list form of the matrix
    """
    aux = []
    for i in list(range(M.shape[1])):
        aux.append(sum(M[:,i]))
    topc = M[:,aux.index(max(aux))]  
    aux[aux.index(max(aux))] = 0
    i = 1
    
    while i<n:
        a = M[:,aux.index(max(aux))] 
        topc = np.column_stack((topc,a))
        aux[aux.index(max(aux))] = 0
        i += 1
    return topc, topc.tolist()
 
    
def get_normalization(M):
    """
    Input: np.array from get_expr_matrix, where rows are genes, columns represent cells
    Action: We perform three types of normalization:
        1. Number of transcripts per 10.000 (we do it for each cell)
        2. log of the above
        3. the above point centered, i.e., zero mean. It's prefered for the PCA 
    Output: three matrices with the normalized gene expressions
    """
    P = []
    for i in range(np.shape(M)[1]):
       P.append(10000*M[:,i]/sum(M[:,i])) 
    P_log = np.copy(np.transpose(np.array(P)))
    for i in range(np.shape(P_log)[0]):
        for j in range(np.shape(P_log)[1]):
            P_log[i,j] = math.log(1+P_log[i,j],10)
    aux = np.copy(P_log)
    for i in range(np.shape(aux)[0]):
        a = aux[i, :]
        #aux[i, :] = (a-np.mean(a))/np.std(a) if np.std(a)>0 else a
        aux[i, :] = a-np.mean(a)
    return np.transpose(np.array(P)), P_log, aux

def get_submatrix(lista, M, G): 
    out = np.copy(M)
    for i in range(np.shape(M)[0]):
       out[i, :] = M[i, :] if G[i] in lista else 0
    return out

def save_gene_list(l, name):
    """    l is the list of genes, name is the file name   """
    f = open(name,'w')
    f.write('\n') # LINE ADDED BY STEFANO
    for i in l:
        f.write(i + '\n')
    f.close()

def save_array(X, name, delim=';'):
    #A = np.array(X)
    A = X
    csv.register_dialect(
        'mydialect',
        delimiter = delim,
        quotechar = '"',
        doublequote = True,
        skipinitialspace = True,
        lineterminator = '\r\n',
        quoting = csv.QUOTE_MINIMAL)
    with open(name, 'w') as mycsvfile:
        thedatawriter = csv.writer(mycsvfile, dialect='mydialect')
        for row in A:
            thedatawriter.writerow(row)  
            
def save_from_list(X, name):
    with open(name, "w") as thefile:
        for item in X:
            thefile.write("%s\n" % item)            
    
"""      
St_genes, St_expr, St_expr_np = get_from_txt('stimulated.txt')
Ct_genes, Ct_expr, Ct_expr_np = get_from_txt('control.txt')
Ct = get_top_cells_matrix(Ct_expr_np,250)[0]
St = get_top_cells_matrix(St_expr_np,250)[0]
Ct_norm = get_normalization(Ct)[0]
St_norm = get_normalization(St)[0]
"""

"""
St_genes, St_expr, St_expr_np = get_from_txt('ATP_1.txt')
Ct_genes, Ct_expr, Ct_expr_np = get_from_txt('CTR_1.txt')
Ct = get_top_cells_matrix(Ct_expr_np,500)[0]
St = get_top_cells_matrix(St_expr_np,500)[0]
Ct_norm = get_normalization(Ct)[0]
St_norm = get_normalization(St)[0]
Ct_log = get_normalization(Ct)[1]
St_log = get_normalization(St)[1]
Ct_sca = get_normalization(Ct)[2]
St_sca = get_normalization(St)[2]
"""



