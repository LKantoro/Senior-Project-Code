from cdlib.classes.node_clustering import NodeClustering
import community.community_louvain as community_louvain
from cdlib import algorithms, viz
import cdlib as cd
from collections import defaultdict

import networkx as nx
import cdlib.algorithms as algo
import cdlib.viz as viz
import matplotlib.pyplot as plt

import numpy as np
import community

import pandas as pd

def louvain(
    g_original: object,
    weight: str = "weight",
    resolution: float = 1.0,
    randomize: int = None) -> NodeClustering:
    g = g_original
    coms = community_louvain.best_partition(
        g, weight=weight, resolution=resolution, randomize=randomize
    )

    # Reshaping the results
    coms_to_node = defaultdict(list)
    for n, c in coms.items():
        coms_to_node[c].append(n)

    coms_louvain = [list(c) for c in coms_to_node.values()]
    return cd.NodeClustering(
        coms_louvain,
        g_original,
        "Louvain",
        method_parameters={
            "weight": weight,
            "resolution": resolution,
            "randomize": randomize,
        },
    )

data = pd.read_excel("Full Pancreatic Cancer Cell Model.xlsx", sheet_name=0)
data.fillna(0, inplace=True)
#print(data)

MatrixA = data
MatrixA = np.array(MatrixA)

MatrixSYM = MatrixA + MatrixA.transpose() - 2*np.diagflat([MatrixA.diagonal()])
# For a symmetric matrix where multiple connections between nodes are added together


MatrixSYM.tolist()

G = nx.from_numpy_array(MatrixSYM)

nx.draw(G, with_labels=True)
plt.show() # symmetric, undirected graph with no self loops

# Can arrange the nodes in a circle
nx.draw_circular(G, with_labels=True)
plt.show() # same thing, but nodes are in circle for easier visualization


node_cluster = louvain(G)

# We extract the communities that produced maximal modularity
com = node_cluster.communities



viz.plot_network_clusters(G, node_cluster, plot_labels=True)


mod = node_cluster.newman_girvan_modularity()
plt.title('Modularity = ' + str(round(mod.score, 5)))
plt.show()
nx.draw(G, with_labels=True)

plt.show()

