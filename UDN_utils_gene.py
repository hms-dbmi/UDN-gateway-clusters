"""TODO description
    TODO import all libraries
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
def get_gene_data(patient_phen, filename, var_or_gene):
    """Retrieve genetic data from a text file (formatted from JSON file)
    Parameters: patient_phen (dict):  dictionary with patients as keys, with values being dictionaries with 
                                        keys ("pos","neg") with a list of the positive and negative 
                                        phenotypes presented by each patient 
                filename (str): name of the text file with the genetic information
                var_or_gene (str): "Var" if variants of "Gen" if genes
    Returns: genomic_data (dict): dictionary with UDN ID as key and list of dictionaries as value, 
                            each dictionary containing information about genes or variants
    """
    genomic_data={}
    with open(filename,"r") as pg:
        lines=pg.readlines()
        for line in lines:
            if line.split("<")[0]=="ID":
                pid=line.split(" ")[3].split("\n")[0]
                genomic_data[pid]=[]
            elif line.split("<")[0]==var_or_gene:
                var=int(line.split(" ")[1].split("\n")[0])
                genomic_data[pid].append({})
            else:
                if not(len(line.split(" "))==1):
                    genomic_data[pid][var][line.split(" ")[0].split("\n")[0]]=line.split(" ")[1].split("\n")[0]
    for patient in genomic_data:
        if not(patient in list(patient_phen.keys())):
            del genomic_data[patient]
    return genomic_data

def get_dist_genomic(genomic_data,var_or_gene):
    """Get the distribution associated to genomic data for its characteristics
    Parameters: genomic_data (dict): with UDN ID as key and list with dictionaries as value, 
                            dict contaning characteristics of the considered genomic data
                var_or_gene (str): "Var" if variants, "Gen" otherwise
    Returns: gene_effects (collec.Counter): distribution of characteristics for selected genomic data
    """
    gene_list=[]
    for patient in genomic_data:
        for i in range(len(genomic_data[patient])):
            if var_or_gene=="Var":
                if "effect" in list(genomic_data[patient][i].keys()) and "gene" in list(genomic_data[patient][i].keys()):
                    gene_list.append([genomic_data[patient][i]["gene"],genomic_data[patient][i]["effect"]])
                else:
                    gene_list.append([genomic_data[patient][i]["gene"],"NA"])
            elif var_or_gene=="Gen":
                if "status" in list(genomic_data[patient][i].keys()) and "gene" in list(genomic_data[patient][i].keys()):
                    gene_list.append([genomic_data[patient][i]["gene"],genomic_data[patient][i]["status"]])
                else:
                    gene_list.append([genomic_data[patient][i]["gene"],"NA"])  
            else:
                print("var_or_gene must be Var or Gen")
    gene_effects=collec.Counter(np.array(gene_list)[:,1])
    return gene_effects
