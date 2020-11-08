#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for SVC

"""

__author__ = "Josephine Yates"
__email__ = "josephine.yates@yahoo.fr"

from UDN_utils import *

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
import pandas as pd

import PicSureHpdsLib
import PicSureClient

from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split

def get_HPO_ids(phenotype_names,mapping_HPO,syn_mapping):
    ids = []
    for term in phenotype_names:
        if term in mapping_HPO:
            ids.append(mapping_HPO[term]["id"])
        else:
            ids.append(syn_mapping[term])
    return ids

def get_diseases_binary(mat_phen_pediatric,mapping_HPO,syn_mapping,logger):
    # get the index of unique phenotypes in the phenotype Dataframe
    header_phen = phenotypes.columns[1:]
    mat_phen_ind=[]
    uniquep=[]
    for i,phen in enumerate(header_phen):
        if not(phen.split("\\")[-2] in uniquep):
            mat_phen_ind.append(i)
            uniquep.append(phen.split("\\")[-2])
    len(mat_phen_ind)
    mat_phen_pediatric=phenotypes.iloc[:,mat_phen_ind]
    mat_phen_pediatric=mat_phen_pediatric.loc[pediatric_patients]
    mat_phen_pediatric=mat_phen_pediatric.replace(to_replace={"Positive": 1, "Negative": 0, np.nan: 0})

    HPO_rett = pd.read_excel("HPO_diseases/HPO_Rett.xlsx", sheet_name=None)["associations"]
    HPO_hurler = pd.read_excel("HPO_diseases/HPO_Hurler.xlsx", sheet_name=None)["associations"]
    HPO_FSHD = pd.read_excel("HPO_diseases/HPO_FSHD.xlsx", sheet_name=None)["associations"]
    HPO_FD = pd.read_excel("HPO_diseases/HPO_familial_dysautonomia.xlsx", sheet_name=None)["associations"]
    HPO_HGPS = pd.read_excel("HPO_diseases/HPO_HGPS.xlsx", sheet_name=None)["associations"]
    HPO_SA = pd.read_excel("HPO_diseases/HPO_spinocerebellar_ataxia_chilhood_onset.xlsx", sheet_name=None)["associations"]
    HPO_MM = pd.read_excel("HPO_diseases/HPO_myopathy_myofibrillar.xlsx", sheet_name=None)["associations"]
    HPO_EE = pd.read_excel("HPO_diseases/HPO_EE_childhood_onset.xlsx", sheet_name=None)["associations"]
    phenotype_names = [mat_phen_pediatric.columns[i].split("\\")[-2] for i in range(len(mat_phen_pediatric.columns))]
    hpo_ids = get_HPO_ids(phenotype_names,mapping_HPO,syn_mapping)
    HPO_dis = {"rett": HPO_rett, "hurler": HPO_hurler,"FSHD": HPO_FSHD,"FD": HPO_FD,"HGPS": HPO_HGPS, "SA": HPO_SA, "MM": HPO_MM, "EE": HPO_EE}
    diseases_binary_HPO = pd.DataFrame(np.zeros((len(list(HPO_dis.keys())),len(phenotype_names))),index=list(HPO_dis.keys()),columns=phenotype_names)
    for dis in diseases_binary_HPO.index:
        for ind in np.where(np.isin(hpo_ids,HPO_dis[dis]["HPO_TERM_ID"].values)):
            diseases_binary_HPO.loc[dis][ind]=1
    return diseases_binary_HPO, mat_phen_pediatric

def SVC_pred(ind_clusters_ped,clusters_ped,mat_phen_pediatric,diseases_binary_HPO,consensus_clustering_labels_ped,logger):
    patients_ped_to_keep = []
    for cl in ind_clusters_ped:
        patients_ped_to_keep+=clusters_ped[cl]
    ind_patients_to_keep = [i for i in range(len(list(mat_phen_pediatric.index))) if list(mat_phen_pediatric.index)[i] in patients_ped_to_keep ]
    cv_score = []
    for fold in range(10):
        X_train, X_test, y_train, y_test = train_test_split(mat_phen_pediatric.iloc[ind_patients_to_keep],consensus_clustering_labels_ped[ind_patients_to_keep],test_size=0.2)
        jac_sim_train = 1 - pairwise_distances(X_train, metric = "jaccard")
        jac_sim_test = 1 - pairwise_distances(X_test, X_train, metric = "jaccard")
        svc = SVC(kernel="precomputed",class_weight="balanced",C=4)
        svc.fit(jac_sim_train,y_train)
        y_pred = svc.predict(jac_sim_test)
        cv_score.append(accuracy_score(y_test,y_pred))
        ## CV 10 fold 0.90
    logger.info("CV 10-fold mean accuracy : {}".format(np.mean(cv_score)))
    jac_sim_train = 1 - pairwise_distances(mat_phen_pediatric.iloc[ind_patients_to_keep], metric = "jaccard")
    dis_sim_ped = 1 - pairwise_distances(diseases_binary_HPO, mat_phen_pediatric.iloc[ind_patients_to_keep], metric = "jaccard")
    svc = SVC(kernel="precomputed",class_weight="balanced",C=4)
    svc.fit(jac_sim_train,consensus_clustering_labels_ped[ind_patients_to_keep])
    y_pred = svc.predict(dis_sim_ped)
    logger.info("The predictions for standard patients of diseases are {}".format(y_pred))
    return svc,y_pred

