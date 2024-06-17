import igraph as ig
import leidenalg as la
import matplotlib.pyplot as plt
import networkx as nx
import cdlib.algorithms as algo
import numpy as np
import community
import cdlib.viz as viz
import pandas as pd



#Reading in data
data = pd.read_excel("Full Pancreatic Cancer Cell Model.xlsx", sheet_name=0)
data.fillna(0, inplace=True)


MatrixA = data
MatrixA = np.matrix(MatrixA)
MatrixA = np.array(MatrixA)
MatrixA.tolist()
MatrixA = pd.DataFrame(MatrixA)


data = MatrixA.values > 0
graph = ig.Graph.Adjacency((data).tolist(), mode = "directed")


graph.es['weight'] = data[data.nonzero()]

partitionlist = []
modularitylist = []
for i in range(1,1000):

    partition = la.find_partition(graph, la.ModularityVertexPartition, seed= i) # takes in a pandas data frame
    cancermodularity = graph.modularity(partition, weights = graph.es['weight'])
    partitionlist.append(partition)
    modularitylist.append(cancermodularity)


from collections import Counter
Counter(modularitylist) #maximum modularity partition is

#cancermodularity = graph.modularity(partition)


# visualization
import itertools
r1=0
r2=68
#using the chain function from the itertools module
numbers = list(itertools.chain(range(r1, r2+1)))


ig.plot(partitionlist[8], "directedleidenpancreaticmodel.png", vertex_size = 20, vertex_label=numbers)
ig.plot(graph,"banana.png", vertex_size = 20, vertex_label = numbers)
