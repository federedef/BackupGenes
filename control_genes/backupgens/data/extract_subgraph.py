#!/usr/bin/env python
# coding: utf-8
import argparse
import networkx as nx
import numpy as np

def get_node_labels(file):
        f = open(file, 'r')
        nodes = []
        for line in f:
                node = line.strip("\n")
                nodes.append(node)
        f.close()
        return nodes

def Adjmatrix2Net(Matrix, rowIds, colIds):
        relations = []
        G = nx.Graph()
        for rowPos, rowId in enumerate(rowIds):
                for colPos, colId in enumerate(rowIds):
                        associationValue = Matrix[rowPos, colPos]
                        if associationValue > 0: G.add_edge(rowId, colId, weight=associationValue)
        return G

def get_seedgroups(file):
        f = open(file, 'r')
        seedgroups = {}
        for line in f:
                fields = line.strip("\n").split("\t")
                seed_name = fields[0]
                seeds = fields[1].split(",")
                seedgroups[seed_name] = seeds
        f.close()
        return seedgroups

if __name__=="__main__":
        parser = argparse.ArgumentParser(description="Add the adcacency matrix and the output name for the embbeded matrix")
        parser.add_argument("-i", "--input", dest="input",
                help="File in numpy format to use as matrix: kernel matrix or adjacency matrix")
        parser.add_argument("-n","--nodes",dest="nodes",required=False,
                help="File with the names of each node")
        parser.add_argument("-s", dest="seedgroups",help="a list with every subgraph needed")
        parser.add_argument("-o", "--output", required=False,dest="output", default="clusters.txt",
                help="Output file")
        options = parser.parse_args()

A=np.load(options.input)

if options.nodes is not None:
        nodes = get_node_labels(options.nodes)
else:
        nodes = list(range(0, A.shape[0] + 1))

graph = Adjmatrix2Net(A, nodes, nodes)

seedgroups = get_seedgroups(options.seedgroups)
for seed_name,seeds in seedgroups.items():
        subgraph = graph.subgraph(seeds)
        nx.write_weighted_edgelist(subgraph, seed_name + "_subgraph", delimiter="\t")