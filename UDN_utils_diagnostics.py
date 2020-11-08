#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for analysis of diagnostics in clusters   

"""

__author__ = "Josephine Yates"
__email__ = "josephine.yates@yahoo.fr"

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
from scipy.stats import mannwhitneyu, chisquare
from sklearn.metrics.pairwise import pairwise_distances
from docx import Document
from docx.shared import Inches
import ast
import logging
from collections import Counter

import PicSureHpdsLib
import PicSureClient

def get_diagnoses(clusters_ad,diagnostics,logger):
    diseases_ad = {cluster: [] for cluster in clusters_ad}
    for cluster in clusters_ad:
        for pat in clusters_ad[cluster]:
            diag = list(diagnostics.loc[pat].values)
            if type(diag)==float:
                continue
            else:
                diseases_ad[cluster]+=diag
    return diseases_ad

def log_diagnoses(diseases_ad,logger):
    for cl in diseases_ad:
        n_dis = 0
        for dis in diseases_ad[cl]:
            if not(type(dis)==float) and not(dis=="affected") and not(dis==""):
                n_dis+=1
        logger.info("Cluster {} has {} diagnosed diseases".format(cl,n_dis))
    
    for cl in diseases_ad:
        logger.info("Cluster ",cl)
        cter= Counter(diseases_ad[cl])
        for dis in np.unique(diseases_ad[cl]):
            logger.info("{}, {}".format(dis,cter[dis]))