"""TODO write description
    TODO comment all functions
"""
from UDN_utils import *
from UDN_utils_cluster_analysis import *
from UDN_utils_clustering import *
from UDN_utils_disease_enrichment import *
from UDN_utils_gene import *
from UDN_utils_HPO_analysis import *
from UDN_utils_parental_age import *
from UDN_utils_primary_symptoms import *

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

def get_link_between_PS_HPO(patient_phen,primary_symptoms,list_phenotypes_unique):
    """Returns the link count of occurrence of a certain HPO large group for patients with a certain primary symptom
    Parameters : patient_phen (dict):  dictionary with patients as keys, with values being dictionaries with 
                                        keys ("pos","neg") with a list of the positive and negative 
                                        phenotypes presented by each patient 
                 primary_symptoms (pandas.core.DataFrame): dataframe with UDN IDs as index, 
                                        and list of primary symptoms reported 
                 list_phenotypes_unique (dict): dictionary of link between phenotypes and the large groups 
                                        they are linked to in the HPO hierarchy

    Returns : link_PS_HPO (dict): with keys ("pos","neg") that contain a dictionary with the 
                    primary symptoms as keys and a dictionary 
                    with the count for every large group of HPO hierarchy of occurrences as value
    """
    link_PS_HPO={"pos": {}, "neg": {}}
    for patient in patient_phen:
        ps=list(primary_symptoms.loc[patient])[1]
        if not(ps in link_PS_HPO["pos"]):
            link_PS_HPO["pos"][ps]={}
        if not(ps in link_PS_HPO["neg"]):
            link_PS_HPO["neg"][ps]={}
        for phen in patient_phen[patient]["pos"]:
            for lg in list_phenotypes_unique[phen]:
                if lg in link_PS_HPO["pos"][ps]:
                    link_PS_HPO["pos"][ps][lg]+=1
                else:
                    link_PS_HPO["pos"][ps][lg]=1
        for phen in patient_phen[patient]["neg"]:
            for lg in list_phenotypes_unique[phen]:
                if lg in link_PS_HPO["neg"][ps]:
                    link_PS_HPO["neg"][ps][lg]+=1
                else:
                    link_PS_HPO["neg"][ps][lg]=1
    return link_PS_HPO

