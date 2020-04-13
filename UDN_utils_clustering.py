"""TODO description & libraries
"""
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import scipy.stats as st
import scipy.sparse as sp
from scipy.stats import fisher_exact
import networkx as nx
from community import community_louvain
from scipy.stats import kruskal
import seaborn as sns
import collections as collec
import os
import xml.etree.ElementTree as ET
import operator
import pandas
import csv
from scipy.stats import mannwhitneyu
from sklearn.metrics.pairwise import pairwise_distances
from docx import Document
from docx.shared import Inches
import ast

import PicSureHpdsLib
import PicSureClient

def get_unique_phenotype(phenotypes):
    """get the list of unique phenotypes in the UDN database
    Parameters: phenotypes (pandas.core.DataFrame): pandas dataframe with the phenotypes

    Returns: mat_phen_ind (list): list of unique phenotypes
    """
    mat_phen_ind, uniquep = [],[]
    header_phen=list(phenotypes)[1:]
    for i,phen in enumerate(header_phen):
        if not(phen.split("\\")[-2] in uniquep):
            mat_phen_ind.append(i)
            uniquep.append(phen.split("\\")[-2])
    return mat_phen_ind

def graph_of_patients_js(UDN_IDs,sim_matrix):
    """Constructs the graph of UDN patients using the similarity matrix computed: nodes are patients, 
                edges between patient i and j is proportional to the similarity between these two patients

    Parameters: UDN_IDs (list): list of UDN IDs of patients to consider
                sim_matrix (np.ndarray): array, similarity matrix of pairwise similarity between each patient

    Returns : G (networkx.Graph): networkx graph of UDN patients 
              pos (np.ndarray): positions of nodes 
    """
    G= nx.Graph()
    elist=[]
    for i in range(sim_matrix.shape[0]):
        G.add_node(UDN_IDs[i])
        for j in range(i,sim_matrix.shape[1]):
            elist.append((UDN_IDs[i],UDN_IDs[j],sim_matrix[i,j]))
    G.add_weighted_edges_from(elist)
    pos=nx.spring_layout(G,dim=2)
    return G,pos

def compute_clusters_community(graph,resolution,logger):
    """Compute the clusters in a graph using Louvain's community detection method
    Parameters : graph (networkx.Graph): graph of UDN patients computed using the pairwise 
                                similarity between patients
                resolution (float): resolution for the Louvain method
                logger (Logger): logger

    Returns: clusters (dict): dictionary with the cluster number as key and a list containing all 
                                the patients in the cluster as value
    """
    partition = community_louvain.best_partition(graph,resolution=resolution)
    logger.info("Partition done")
    clusters={}
    for node in partition.keys():
        if not(partition[node] in clusters.keys()):
            clusters[partition[node]]=[node]
        else:
            clusters[partition[node]].append(node)
    count=0
    for cluster in clusters.keys():
        logger.info("Length of cluster {} : {}".format(cluster,len(clusters[cluster])))
        if len(clusters[cluster])<=1:
            count+=1
    logger.info("Number of clusters with less than 3 patients (outliers) : {}".format(count))
    return clusters

