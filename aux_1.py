#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 23:30:14 2016

@author: tomasz.ignac
"""

#import genes_1
from sklearn.decomposition import NMF
from sklearn.cluster import KMeans

def MI(M,d):
    out = np.empty((2,np.shape(M)[0]))
    for row in range(np.shape(M)[0]):
        g = M[row,:]
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
        out[:,row] = mi, gt
    return out

def save_array(X):
    A = np.array(X)
    c = 0
    csv.register_dialect(
        'mydialect',
        delimiter = ';',
        quotechar = '"',
        doublequote = True,
        skipinitialspace = True,
        lineterminator = '\r\n',
        quoting = csv.QUOTE_MINIMAL)
    with open('mydata.csv', 'w') as mycsvfile:
        thedatawriter = csv.writer(mycsvfile, dialect='mydialect')
        for row in A:
            thedatawriter.writerow(row)
            c = c+1
    return c

def get_prediction_loop(X,Y,tr=0.8,it=1000,d=None):
    if d is not None:
        X[d] = [0]*len(X[d])
        Y[d] = [0]*len(Y[d])
    pair = np.hstack((X, Y))
    b=[]
    a = get_binary(pair)
    B = Compare((a[:,:np.shape(X)[1]], Genes),(a[:,np.shape(X)[1]:np.shape(pair)[1]], Genes)) 
    b.append(B.LR_loop(Genes,tr,it))
    b.append(B.NB_loop(Genes,tr,it))
    b.append(B.RF_loop(Genes,tr,it))
    C = Compare((X,Genes),(Y,Genes))
    b.append(C.LR_loop(Genes,tr,it))
    b.append(C.NB_loop(Genes,tr,it,False))
    b.append(C.NB_loop(Genes,tr,it))
    b.append(C.RF_loop(Genes,tr,it,False))
    b.append(C.RF_loop(Genes,tr,it))
    return b,B,C

def aux_fcol(X,Y,Z=None):
    aux = (X,Y,Z)
    out = []
    i = 2 if Z is None else 3
    for i in range(i):
        a = np.shape(aux[i])[1]*['blue']
        a[:50] = 50*['red']
        a[-50:] = 50*['green']
        out = out + a
    return out

    
def get_DR_bis(X, Y, com, Z=None):
    pair = np.hstack((X, Y, Z)) if Z is not None else np.hstack((X, Y))
    y_color = np.shape(X)[1]*['red'] + np.shape(Y)[1]*['blue']
    y_color = y_color + np.shape(Z)[1]*['green'] if Z is not None else y_color
    #y_color = aux_fcol(X,Y,Z)
    tsne = TSNE(n_components=com, random_state=0)   
    y_tsne = tsne.fit_transform(np.transpose(pair))  
    pca = PCA(n_components=com)
    y_pca = pca.fit_transform(np.transpose(pair))
    if com == 2:
        plt.figure()
        plt.scatter(y_pca[:,0],y_pca[:,1], c=y_color)
        plt.figure()
        plt.scatter(y_tsne[:,0],y_tsne[:,1], c=y_color)
    elif com == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(y_pca[:,0], y_pca[:,1], y_pca[:,2], c=y_color)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(y_tsne[:,0], y_tsne[:,1], y_tsne[:,2], c=y_color)
    return y_pca, y_tsne

    
def get_table_19genes():
    o1, o2, o3 = [], [], []
    for i in range(3):
        o1.append(get_prediction_loop(ERY[i],MYL[i])[0])
        o2.append(get_prediction_loop(ERY[i],COM[i])[0])
        o3.append(get_prediction_loop(MYL[i],COM[i])[0])
    out = o1 + o2 +o3
    return out

    
def aux_f1(X,Y,med):
    if med is True:
        pair = np.hstack((X, Y))
        a = get_binary(pair)
        A = Compare((a[:,:np.shape(X)[1]], Genes),(a[:,np.shape(X)[1]:np.shape(pair)[1]], Genes)) 
    else:
        A = Compare((X,Genes), (Y,Genes))
    out_1 = []
    out_2 = []
    for g in Genes:
        aux = A.MI(g)
        out_1.append(aux[0])
        out_2.append(aux[1]) 
    return out_1, out_2
  
    
def get_table_1gene(med=True):
    out = []
    for i in range(3):
        X = ERY[i]
        Y = MYL[i]
        a, b = aux_f1(X,Y,med)    
        out.append(a)
        out.append(b)
    for i in range(3):
        X = ERY[i]
        Y = COM[i]
        a, b = aux_f1(X,Y,med)    
        out.append(a)
        out.append(b)
    for i in range(3):
        X = MYL[i]
        Y = COM[i]
        a, b = aux_f1(X,Y,med)    
        out.append(a)
        out.append(b)
    return out


def get_Distance(X,Y=None,Z=None,ham=False,plot=False):
    if Y is not None and Z is not None:
        M = np.hstack((X, Y, Z))
    elif Y is not None:
        M = np.hstack((X, Y))
    else:
        M = X
    M = get_binary(M) if ham is not False else M
    M = squareform(pdist(np.transpose(M))) if ham is False else squareform(pdist(np.transpose(M), 'hamming')) 
    if plot is not False:
        plt.pcolor(M)
        plt.colorbar()
    return M
   
def bootstrap(v,t):
    out = np.empty((2,t))
    for i in range(t):
        a = np.random.choice(v,500)
        b = np.random.choice(v,500)
        out[0,i] = np.percentile(a,95)
        out[1,i] = abs(np.percentile(b,95)-np.percentile(a,95))
    return out

def pair_rand(v,v1,t):
    out = np.zeros(t)
    for i in range(t):
        a = np.random.choice(v,np.shape(v)[0],replace=False)
        e = 0
        for j in range(np.shape(a)[0]):
            e = e + v1[j] if a[j] == 1 else  e
        out[i] = 100*e/np.sum(v1) if e>0 else 0
    return out

def get_score(lista, M, G): 
    out = np.zeros(np.shape(M)[1])
    ind = np.zeros(np.shape(M)[0])
    for i in lista:
        try:
            ind[G.index(i)] = 1
        except ValueError:
            pass
    for i in range(np.shape(M)[1]):
        out[i] = sum(ind*M[:,i])
    return out
    
def get_submatrix(lista, M, G): 
    out = np.copy(M)
    for i in range(np.shape(M)[0]):
       out[i, :] = M[i, :] if G[i] in lista else 0
    return out
 #VIM, odpowiada za tranzycje 

#pair = np.copy(St_sca)
#pair[21335,:] = 0 #TAOK1 ST
#pair[10639,:] = 0 #Ct, 
#pair[11317,:] = 0 #St, 
gene_c = Ct_genes.index('VIM')
gene_s = St_genes.index('VIM')
#v = np.hstack((Ct_norm[gene,:],St_norm[gene_s,:]))

#plt.style.use('bmh')
#plt.hist(St_norm[gene_s,:],bins=40, histtype="stepfilled", label='Stimulated', stacked=True)
#plt.hist(Ct_norm[gene,:],bins=11, histtype="step", label='Controls', stacked=True)
#plt.legend()
#plt.title('MT-RNR2')
#plt.xlabel('Number of transcripts per 10 000')
#plt.show()

com = 2
tsne = TSNE(n_components=com, random_state=0)   
y_tsne = tsne.fit_transform(np.transpose(pair))  
pca = PCA(n_components=com)
y_pca = pca.fit_transform(np.transpose(pair))
if com == 2:
    plt.figure()
    plt.scatter(y_pca[:,0],y_pca[:,1], c=St_log[gene_s,:])
    plt.colorbar()
    plt.title('Stimulated - PCA\npoints colored by expression of VIM')
    plt.figure()
    plt.scatter(y_tsne[:,0],y_tsne[:,1], c=St_log[gene_s,:])
    plt.title('Stimulated - tSNE\npoints colored by expression of VIM')
    plt.colorbar()
elif com == 3:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(y_pca[:,0], y_pca[:,1], y_pca[:,2])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(y_tsne[:,0], y_tsne[:,1], y_tsne[:,2])

def get_top_genes(M,names,n):
    out = []
    aux = np.empty(np.shape(M)[0])
    for i in range(np.shape(M)[0]):
        aux[i] = np.sum(M[i,:])
    #aux = np.argsort(aux)
    for i in range(n):
        p = np.argmax(aux)
        out.append(names[p])
        aux[p] = 0
    return out

"""
def get_corr(M,p):
    out = np.empty(np.shape(M)[0])
    for i in range(np.shape(M)[0]):
        a = np.vstack((M[p,:],M[i,:]))
        aux = np.corrcoef(a)[0,1]
        out[i] = 0 if np.isnan(aux) else aux
    return out
a = get_corr(St_norm, St_genes.index('MT-RNR2'))  
def get_reduced(M,a,per):
    out = np.copy(M)
    aux = np.percentile(a,per)
    for i in range(np.shape(a)[0]):
        out[i,:] = 0 if a[i]>aux else out[i,:]
    return out
  


    
    

#kmeans = KMeans(n_clusters=5, random_state=0).fit(y_tsne)    
#aux=kmeans.labels_
#plt.scatter(y_tsne[:,0],y_tsne[:,1], c=aux) 
#plt.colorbar()
#plt.figure()
#plt.scatter(y_pca[:,0],y_pca[:,1], c=aux)   
#plt.colorbar()


aux = np.zeros(500)
for i in range(500):
    #if y_tsne[i,0] < -10 and y_tsne[i,0] > -30 and y_tsne[i,1]< 0:
    if y_pca[i,0] < 1 and y_pca[i,0] > -0.5:
        aux[i] = 1
        
#g = Ct_genes.index('VIM')
#kexpr = np.zeros(5)
#for i in range(np.shape(aux)[0]):
#    kexpr[aux[i]] = kexpr[aux[i]] + Ct_norm[g,i]

expr_2 = []
for i in St_genes:    
    g = St_genes.index(i)
    e = 0
    for j in range(np.shape(aux)[0]):
        e = e + St_log[g,j] if aux[j] == 1 else  e
    e = 100*e/np.sum(St_log[g,:]) if e>0 else 0
    expr_2.append([e, min(np.sum(St_norm[g,:]),1000)])
del(e)
expr_2=np.array(expr_2)   
ind = []
for i in range(np.shape(expr_2)[0]):
    if expr_2[i,0]>23 and expr_2[i,1]>400:
        ind.append(St_genes[i])



       #lista_st.append([St_genes[i], mi_st_tsne[0,i]])
        
#plt.pcolor(aux)
#plt.colorbar()
    
#fig = plt.figure()
#fig.add_subplot(2, 2, 1)
#plt.plot(aux_1[0,:],marker='^',color='red',ls='dotted')
#fig.add_subplot(2, 2, 2)
#plt.plot(aux_1[1,:],marker='o',color='blue',ls='--')

#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#plt.plot(aux_1[0,:],marker='^',color='red',ls='dotted')
#plt.plot(aux_1[1,:],marker='o',color='blue',ls='--')


#plt.hist(X[0],bins=25)
#plt.xlabel('Expr')

#plt.show()

#model = NMF(n_components=2, init='random', random_state=0)
#M=np.copy(Ct_norm)
#M[Ct_genes.index('MT-RNR2'),:] = 0
#M=get_reduced(St_norm,a,99)
#g = St_genes.index('MT-RNR2' )
#y = model.fit_transform(np.transpose(M)) 
#y = np.transpose(y)
#plt.figure()
#plt.scatter(y[0,:],y[1,:],  c=St_norm[g,:])
#plt.colorbar()


#plt.figure()
#plt.scatter(y_pca[:,0],y_pca[:,1], c=SC[2780,:])
#plt.colorbar()
#plt.figure()
#plt.scatter(y_tsne[:,0],y_tsne[:,1], c=SC[2780,:])
#plt.colorbar()
"""