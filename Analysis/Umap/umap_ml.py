<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 15:02:21 2019

@author: Parthiv_Lab
"""

import sys
import pandas
import numpy
import matplotlib.pyplot as plt
import umap
import pickle

import statsmodels.api as sm
from sklearn.cluster import KMeans
from pynndescent import NNDescent
import community
import networkx as nx


## Import and embedding

out = pandas.read_csv('mixout.csv', header=None).values.squeeze()
mtx = pandas.read_csv('mixall.csv', header=None).values.squeeze()


numpy.random.seed(0)
mdl = KMeans(init='k-means++', n_clusters=30, n_init=10)
kmeansfit = mdl.fit_predict(mtx[:,3:15])

numpy.random.seed(0)
reducer = umap.UMAP(n_neighbors=12,n_components=2,metric='chebyshev', min_dist=.1)
embedding = reducer.fit_transform(mtx[:,3:15])
plt.scatter(embedding[:, 0], embedding[:, 1],3,c=kmeansfit,cmap='hsv')


#numpy.savetxt("embeddinumap.csv", embedding, delimiter=",")
embedding = pandas.read_csv('embeddingumap.csv', header=None).values.squeeze()


## Ajacency matrix

NNindex = NNDescent(embedding)
adjacency = NNindex.query(embedding, k=15)

adj_matrix = []
weighted_matrix = [] 
for i in range(0,len(embedding[:,1])):
    temp_adj = numpy.zeros(len(embedding)).tolist()
    temp_weight = numpy.zeros(len(embedding)).tolist()
    for j in range(1,len(adjacency[0][i])):
        temp_adj[adjacency[0][i][j]]=1
        temp_weight[adjacency[0][i][j]]=adjacency[1][i][j]
    adj_matrix.append(temp_adj)        
    weighted_matrix.append(temp_weight)

del temp_adj
del temp_weight
weighted_matrix = numpy.array(weighted_matrix)
adj_matrix = numpy.array(adj_matrix)


## Create graph
G = nx.from_numpy_matrix(weighted_matrix)

#first compute the best partition
partition = community.best_partition(G)
#dendo1 = community.generate_dendrogram(G)
#testcomm = community.induced_graph(dendo1[0], G)

coloringlist =  list(partition.values())
plt.scatter(embedding[:,0], embedding[:,1], c=coloringlist, cmap = 'hsv',s=8,alpha = .2)




coloringlist = numpy.array(coloringlist)
## Define cell groups and find centroid
groupindex = []
groupindex_trunc = []
meanembedding = []
for i in range(0,max(coloringlist)+1):
    classestrunc_temp = []
    classes_temp = numpy.where(coloringlist==i)[0]
    for j in classes_temp:
        if j <= len(out[:,1]):
            classestrunc_temp.append(j)
            
    meanembedding.append([numpy.mean(embedding[classes_temp,0]),numpy.mean(embedding[classes_temp,1])])
    groupindex.append(classes_temp)
    groupindex_trunc.append(classestrunc_temp)
    
meanembedding = numpy.array(meanembedding)



plt.scatter(embedding[:,0], embedding[:,1], c=coloringlist, cmap = 'hsv',s=8,alpha = .2) #,c=mtx[:,14],s=4,alpha = .1,cmap = 'inferno'
for i in range(0,max(coloringlist)+1):
    plt.annotate(range(0,max(coloringlist)+1)[i], (meanembedding[i,0], meanembedding[i,1]))






meanlist = []
for i in groupindex_trunc:
    templist = out1[i,:]    
    gmeanact = []
    for j in range(1,5):
        #gmeanact.append(numpy.mean(numpy.log(templist[numpy.where(templist[:,0]==j),1])))
        gmeanact.append(numpy.mean(templist[numpy.where(templist[:,0]==j),3]))
    meanlist.append(gmeanact)

for i in meanlist:
    plt.plot(range(0,4),i) #,c=mtx[:,14],s=4,alpha = .1,cmap = 'inferno'

meanlist = numpy.array(meanlist)
=======
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 15:02:21 2019

@author: Parthiv_Lab
"""

import sys
import pandas
import numpy
import matplotlib.pyplot as plt
import umap
import pickle

import statsmodels.api as sm
from sklearn.cluster import KMeans
from pynndescent import NNDescent
import community
import networkx as nx


## Import and embedding

out = pandas.read_csv('mixout.csv', header=None).values.squeeze()
mtx = pandas.read_csv('mixall.csv', header=None).values.squeeze()


numpy.random.seed(0)
mdl = KMeans(init='k-means++', n_clusters=30, n_init=10)
kmeansfit = mdl.fit_predict(mtx[:,3:15])

numpy.random.seed(0)
reducer = umap.UMAP(n_neighbors=12,n_components=2,metric='chebyshev', min_dist=.1)
embedding = reducer.fit_transform(mtx[:,3:15])
plt.scatter(embedding[:, 0], embedding[:, 1],3,c=kmeansfit,cmap='hsv')


#numpy.savetxt("embeddinumap.csv", embedding, delimiter=",")
embedding = pandas.read_csv('embeddingumap.csv', header=None).values.squeeze()


## Ajacency matrix

NNindex = NNDescent(embedding)
adjacency = NNindex.query(embedding, k=15)

adj_matrix = []
weighted_matrix = [] 
for i in range(0,len(embedding[:,1])):
    temp_adj = numpy.zeros(len(embedding)).tolist()
    temp_weight = numpy.zeros(len(embedding)).tolist()
    for j in range(1,len(adjacency[0][i])):
        temp_adj[adjacency[0][i][j]]=1
        temp_weight[adjacency[0][i][j]]=adjacency[1][i][j]
    adj_matrix.append(temp_adj)        
    weighted_matrix.append(temp_weight)

del temp_adj
del temp_weight
weighted_matrix = numpy.array(weighted_matrix)
adj_matrix = numpy.array(adj_matrix)


## Create graph
G = nx.from_numpy_matrix(weighted_matrix)

#first compute the best partition
partition = community.best_partition(G)
#dendo1 = community.generate_dendrogram(G)
#testcomm = community.induced_graph(dendo1[0], G)

coloringlist =  list(partition.values())
plt.scatter(embedding[:,0], embedding[:,1], c=coloringlist, cmap = 'hsv',s=8,alpha = .2)




coloringlist = numpy.array(coloringlist)
## Define cell groups and find centroid
groupindex = []
groupindex_trunc = []
meanembedding = []
for i in range(0,max(coloringlist)+1):
    classestrunc_temp = []
    classes_temp = numpy.where(coloringlist==i)[0]
    for j in classes_temp:
        if j <= len(out[:,1]):
            classestrunc_temp.append(j)
            
    meanembedding.append([numpy.mean(embedding[classes_temp,0]),numpy.mean(embedding[classes_temp,1])])
    groupindex.append(classes_temp)
    groupindex_trunc.append(classestrunc_temp)
    
meanembedding = numpy.array(meanembedding)



plt.scatter(embedding[:,0], embedding[:,1], c=coloringlist, cmap = 'hsv',s=8,alpha = .2) #,c=mtx[:,14],s=4,alpha = .1,cmap = 'inferno'
for i in range(0,max(coloringlist)+1):
    plt.annotate(range(0,max(coloringlist)+1)[i], (meanembedding[i,0], meanembedding[i,1]))






meanlist = []
for i in groupindex_trunc:
    templist = out1[i,:]    
    gmeanact = []
    for j in range(1,5):
        #gmeanact.append(numpy.mean(numpy.log(templist[numpy.where(templist[:,0]==j),1])))
        gmeanact.append(numpy.mean(templist[numpy.where(templist[:,0]==j),3]))
    meanlist.append(gmeanact)

for i in meanlist:
    plt.plot(range(0,4),i) #,c=mtx[:,14],s=4,alpha = .1,cmap = 'inferno'

meanlist = numpy.array(meanlist)
>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
