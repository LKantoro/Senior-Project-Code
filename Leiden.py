import igraph as ig
import leidenalg as la
import numpy as np
import pandas as pd

#Reading in data
data = pd.read_excel("Full Pancreatic Cancer Cell Model.xlsx", sheet_name=0)
data.fillna(0, inplace=True)


MatrixA = data
MatrixA = np.matrix(MatrixA)
MatrixA = np.array(MatrixA)

# For a symmetric matrix where multiple connections between nodes are added together
MatrixSYM = MatrixA + MatrixA.transpose() - 2*np.diagflat([MatrixA.diagonal()])
MatrixSYM = MatrixSYM
MatrixSYM.tolist()

#Undirected

MatrixSYM = pd.DataFrame(MatrixSYM)

symmdata = MatrixSYM.values > 0
symmgraph = ig.Graph.Adjacency((symmdata).tolist(), mode = "undirected")

symmgraph.es['weight'] = symmdata[symmdata.nonzero()]



partition = la.find_partition(symmgraph, la.ModularityVertexPartition, seed= 35) # takes in a pandas data frame
cancermodularity = symmgraph.modularity(partition, weights = symmgraph.es['weight'])


partitionlist = []
modularitylist = []
for i in range(1,1000):

    partition = la.find_partition(symmgraph, la.ModularityVertexPartition, seed= i) # takes in a pandas data frame
    cancermodularity = symmgraph.modularity(partition, weights = symmgraph.es['weight'])
    partitionlist.append(partition)
    modularitylist.append(cancermodularity)


from collections import Counter
Counter(modularitylist) #maximum modularity partition is 0.6402739

# visualization
import itertools
r1=0
r2=68
#using the chain function from the itertools module
numbers = list(itertools.chain(range(r1, r2+1)))


ig.plot(partitionlist[8], "pancreaticmodel.png", vertex_size = 20, vertex_label=numbers)
