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
from sklearn.metrics import accuracy_score, hamming_loss, balanced_accuracy_score, f1_score
from sklearn.model_selection import train_test_split, GridSearchCV, KFold

import pickle

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
    mat_to_keep, consensus_to_keep = mat_phen_pediatric.iloc[ind_patients_to_keep],consensus_clustering_labels_ped[ind_patients_to_keep]
    
    X_train, X_test, y_train, y_test = train_test_split(mat_phen_pediatric.iloc[ind_patients_to_keep],consensus_clustering_labels_ped[ind_patients_to_keep],test_size=0.2)
    
    cv_score = {}
    for C in np.linspace(0.01,1,5):
        cv_score[C] = []
        for fold_index, test_index in kf.split(X_train):

            X_fold, X_val, y_fold, y_val = X_train.values[fold_index],X_train.values[test_index], \
                                            y_train[fold_index], \
                                            y_train[test_index]
            jac_sim_fold = 1 - pairwise_distances(X_fold, metric = "jaccard")
            jac_sim_val = 1 - pairwise_distances(X_val, X_fold, metric = "jaccard")
            svc = SVC(kernel="precomputed",class_weight="balanced",C=C)
            svc.fit(jac_sim_fold,y_fold)
            y_pred = svc.predict(jac_sim_val)
            cv_score[C].append(accuracy_score(y_val,y_pred))
        logger.info("Mean CV 10-fold score for C = {}: {}".format(C,np.mean(cv_score[C])))
        
    jac_sim_train_small = 1 - pairwise_distances(X_train, metric = "jaccard")
    jac_sim_ped_test_train = 1 - pairwise_distances(X_test, X_train, metric = "jaccard")
    svc = SVC(kernel="precomputed",class_weight="balanced",C=0.75)
    svc.fit(jac_sim_train_small,y_train)
    y_pred = svc.predict(jac_sim_ped_test_train)
    logger.info("Hamming loss jaccard svc: {}".format(hamming_loss(y_test,y_pred)))
    logger.info("Accuracy jaccard svc: {}".format(accuracy_score(y_test,y_pred)))
    logger.info("Balanced Accuracy jaccard svc: {}".format(balanced_accuracy_score(y_test,y_pred)))
    
    jac_sim_train = 1 - pairwise_distances(mat_to_keep, metric = "jaccard")
    dis_sim_ped = 1 - pairwise_distances(diseases_binary_HPO, X_train, metric = "jaccard")
    y_pred_dis = svc.predict(dis_sim_ped)
    logger.info("The predictions for standard patients of diseases are {}".format(y_pred_dis))
    filename = 'svc_model.sav'
    pickle.dump(svc, open(filename, 'wb'))
    return svc,y_pred_dis

def SVC_pred_hpo(ind_clusters_ped,clusters_ped,mat_phen_pediatric,diseases_binary_HPO,consensus_clustering_labels_ped,logger):
    patients_ped_to_keep = []
    for cl in ind_clusters_ped:
        patients_ped_to_keep+=clusters_ped[cl]
    ind_patients_to_keep = [i for i in range(len(list(mat_phen_pediatric.index))) if list(mat_phen_pediatric.index)[i] in patients_ped_to_keep ]
    mat_to_keep, consensus_to_keep = mat_phen_pediatric.iloc[ind_patients_to_keep],consensus_clustering_labels_ped[ind_patients_to_keep]
    
    X_train, X_test, y_train, y_test = train_test_split(mat_phen_pediatric.iloc[ind_patients_to_keep],consensus_clustering_labels_ped[ind_patients_to_keep],test_size=0.2)
    
    cv_score = {}
    parameters_svc = {"kernel":("rbf","poly"), "C":np.linspace(0.75,3,5),"class_weight":["balanced"],"gamma":["scale"]}
    svc_hpo = SVC()
    clf = GridSearchCV(svc_hpo,parameters_svc,scoring="balanced_accuracy",cv=10)
    clf.fit(X_train.values,y_train)
    y_pred_hpo = clf.predict(X_test)
    
    logger.info("Hamming loss hpo svc: {}".format(hamming_loss(y_test,y_pred_hpo)))
    logger.info("Accuracy hpo svc: {}".format(accuracy_score(y_test,y_pred_hpo)))
    logger.info("Balanced Accuracy hpo svc: {}".format(balanced_accuracy_score(y_test,y_pred_hpo)))
    logger.info("F1 macro hpo svc: {}".format(f1_score(y_test,y_pred_hpo,average="macro")))
    logger.info("F1 weighted hpo svc: {}".format(f1_score(y_test,y_pred_hpo,average="weighted")))
    
    y_pred_hpo_dis = clf.predict(diseases_binary_HPO)
    logger.info("The predictions for standard patients of diseases are {}".format(y_pred_hpo_dis))
    filename = 'svc_model_hpo.sav'
    pickle.dump(svc, open(filename, 'wb'))
    return clf.best_estimator_,y_pred_hpo_dis



