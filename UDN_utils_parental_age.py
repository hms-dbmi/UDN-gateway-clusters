"""TODO description
    TODO import all libraries
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

def get_diff_parent_age(natal_history_adult, natal_history_pediatric, natalhist, logger):
    """Logs statistical difference in parental age between adult and pediatric (Mann-Whitney U test)
    Parameters : natal_history_adult (pandas.core.DataFrame): natal history for adult patients
                    natal_history_pediatric (pandas.core.DataFrame): natal history for pediatric patients
                    natalhist (str): string to study maternal or paternal age
                    logger (Logger): logger
    Returns : None
    """
    nhp=np.array(list(natal_history_pediatric[natalhist]))
    nhp=[nhp[i] for i in range(len(nhp)) if not(np.isnan(nhp)[i])]
    nha=np.array(list(natal_history_adult[natalhist]))
    nha=[nha[i] for i in range(len(nha)) if not(np.isnan(nha)[i])]
    logger.info("Difference between maternal age for pediatric and adult onset : {}".format(mannwhitneyu(nhp,nha)))

def distrib_age(parent_age, known_dist,tranches,boundaries,mat_or_pat):
    """Shows the distribution of maternal age compared between UDN and the US in 2009
    Parameters : parent_age (np.ndarray): array of parental age in the UDN database
                 known_dist (list): list, known distribution of parental age for age groups given in splits
                 tranches (list): list of str, age groups that correspond to the known distribution 
                 boundaries (np.ndarray): array of 2-D arrays, with the boundaries in int corresponding to the splits given in tranches
                 mat_or_pat (str): "maternal" or "paternal", for maternal or paternal age
    Returns: distrib (dict): dictionary with age distribution in the UDN 
    Shows a joint plot of UDN distribution and known distribution of maternal age
    """
    count_age={}
    for age in parent_age:
        if age in count_age:
            count_age[age]+=1
        else:
            count_age[age]=1
    distrib_age=[0 for i in range(len(tranches))]
    for age in count_age:
        for i in range(len(boundaries)):
            if age>=boundaries[i][0] and age<=boundaries[i][1]:
                distrib_age[i]+=count_age[age]/len(parent_age)*100
    plt.figure(figsize=(20,15))
    plt.plot(tranches,distrib_age,'b',label="Distribution in UDN")
    plt.plot(tranches,known_dist,'r',label="Distribution in USA in 2009")
    plt.xlabel(mat_or_pat+" age at birth",fontsize=20)
    plt.ylabel("Distribution in UDN vs USA in 2009 (%)",fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.,fontsize=20)
    plt.show()
    return distrib_age
