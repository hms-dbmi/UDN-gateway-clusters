#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for HPO analysis

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

import PicSureHpdsLib
import PicSureClient

def get_patient_phenotypes(phenotypes):
    """Gets the list of unique phenotypes presented by the patients of the UDN 
    Parameters: phenotypes (pandas.core.DataFrame): pandas dataframe with the phenotypes
    
    Returns : patient_phen (dict) : dictionary with patients as keys, with values \
                                    dictionaries with keys ("pos","neg") \
                             with a list of the positive and negative phenotypes presented by each patient
    """
    header_phen=list(phenotypes)

    patient_phen={patient: {"pos": [], "neg": []} for patient in phenotypes.index.values}

    for patient,row in phenotypes.iterrows():
        for i,phen in enumerate(row):
            if phen=="Positive":
                if not header_phen[i].split("\\")[-2] in patient_phen[patient]["pos"]:
                    patient_phen[patient]["pos"].append(header_phen[i].split("\\")[-2])
            elif phen=="Negative":
                if not header_phen[i].split("\\")[-2] in patient_phen[patient]["neg"]:
                    patient_phen[patient]["neg"].append(header_phen[i].split("\\")[-2])

    for patient in patient_phen:
        if len(patient_phen[patient]["pos"])==0 and len(patient_phen[patient]["neg"])==0:
            del patient_phen[patient]

    return patient_phen

def get_patient_eval_date(jsonfile, patient_phen):
    """Gets the list of patients with their evaluation date
    Parameters: jsonfile (str): path to the json file
                patient_phen (dict) : dictionary with patients as keys, with values \
                                    dictionaries with keys ("pos","neg") \
                             with a list of the positive and negative phenotypes presented by each patient
    
    Returns : patient_eval_date (dict): dict with patient as key, evaluation date as value
    """
    patient_eval_date = {}

    with open("patient_eval_date.txt","r") as evaldate:
        lines=evaldate.readlines()
        for line in lines:
            patient_eval_date[line.split(" ")[3]]=line.split(" ")[5].split("\n")[0]

    for pat in list(patient_eval_date.keys()):
        if not (pat in list(patient_phen)):
             del patient_eval_date[pat]

    print("Number of patients with no information on eval date : ",collec.Counter([patient_eval_date[pat] for pat in patient_eval_date])["None"])
    
    return patient_eval_date 

def patient_eval_before_2015(patient_eval_date, patient_phen):
    """Gets an updated dictionary of patient phenotype, with patients before 2015 with no negative \
        values (cf paper for explanation of possible bias)
    Parameters: patient_eval_date (dict): dict with patient as key, evaluation date as value
                patient_phen (dict) : dictionary with patients as keys, with values \
                                    dictionaries with keys ("pos","neg") \
                             with a list of the positive and negative phenotypes presented by each patient
    
    Returns : patient_phen_wo_2015 (dict): patient_phen dict with updated negative phenotypes
    """
    list_pat_before_2015=[]

    for pat in patient_eval_date:
        if patient_eval_date[pat]=="None":
            continue
        if int(patient_eval_date[pat].split("-")[0])<=2015:
            list_pat_before_2015.append(pat)

    patient_phen_wo_2015=patient_phen.copy()

    for pat in list_pat_before_2015:
        patient_phen_wo_2015[pat]["neg"]=[]

    return patient_phen_wo_2015

def get_HPO_count_list(patient_phen, patient_list):
    """Gets a list of HPO counts for all patients in the patient list 
    Parameters: patient_list (np.ndarray): list of patients to consider
                patient_phen (dict) : dictionary with patients as keys, with values \
                                    dictionaries with keys ("pos","neg") \
                             with a list of the positive and negative phenotypes presented by each patient
    
    Returns : HPO_list_pos, HPO_list_neg, HPO_list (np.ndarray): lists of the count of positive, \
        negative and all counts of HPO terms resp.
    """
    HPO_list_pos = np.array([len(list(patient_phen[patient]["pos"].values())) for patient in patient_list])
    HPO_list_neg = np.array([len(list(patient_phen[patient]["neg"].values())) for patient in patient_list])
    HPO_list = HPO_list_pos + HPO_list_neg

    return HPO_list_pos, HPO_list_neg, HPO_list

def get_HPO_terms(patient_phen, patient_list):
    """Gets HPO count per patient
    Parameters: patient_list (np.ndarray): list of patients to consider
                patient_phen (dict) : dictionary with patients as keys, with values \
                                    dictionaries with keys ("pos","neg") \
                             with a list of the positive and negative phenotypes presented by each patient
    
    Returns : HPO_terms (dict): dictionary with patient as key and HPO term count as value
    """
    HPO_terms = {patient: len(list(patient_phen[patient]["pos"].values())) + \
                        len(list(patient_phen[patient]["neg"].values())) for patient in patient_list }

    return HPO_terms

def get_best_phenotypes(list_patients, patient_phen, nb_of_phen, logger):
    """Shows the phenotypes the most represented in the UDN gateway for a given community of patients
        Parameters: list_patients (list): list of patients IDs that should be considered
                    patient_phen (dict): dictionary with patients as keys, 
                                        with values being dictionaries with keys ("pos","neg") 
                                        with a list of the positive and negative phenotypes presented by each patient
                    nb_of_phen (int): number of best phen to represent
                    logger (Logger): logger

        Returns: list_neg_phen (np.ndarray): list of ranked best negative associations
                 list_pos_phen (np.ndarray): list of ranked best positive associations
        Shows the nb_of_phen best pos and neg HPO association with % of representation
    """
    list_neg_phen,list_pos_phen=[],[]
    neg_phen_w_count,pos_phen_w_count={},{}
    for patient in list_patients:
        for phen in patient_phen[patient]["neg"]:
            if not(phen in neg_phen_w_count):
                neg_phen_w_count[phen]=1/len(list_patients)*100
            else:
                neg_phen_w_count[phen]+=1/len(list_patients)*100
        for phen in patient_phen[patient]["pos"]:
            if not(phen in pos_phen_w_count):
                pos_phen_w_count[phen]=1/len(list_patients)*100
            else:
                pos_phen_w_count[phen]+=1/len(list_patients)*100  
    sorted_dict_pos=collec.OrderedDict(sorted(pos_phen_w_count.items(), key=operator.itemgetter(1), reverse=True))
    sorted_dict_neg=collec.OrderedDict(sorted(neg_phen_w_count.items(), key=operator.itemgetter(1), reverse=True))
    logger.info("Most highly ranked positive phenotypes")
    for i,key in enumerate(sorted_dict_pos):
        logger.info(key,sorted_dict_pos[key])
        list_pos_phen.append(sorted_dict_pos[key])
        if i>nb_of_phen:
            break
    logger.info("Most highly ranked negative phenotypes")
    for i,key in enumerate(sorted_dict_neg):
        logger.info(key,sorted_dict_neg[key])
        list_neg_phen.append(sorted_dict_neg[key])
        if i>nb_of_phen:
            break
    return np.array(list_neg_phen),np.array(list_pos_phen)

def show_stats_HPO_counts(HPO_list,HPO_list_pos,HPO_list_neg, logger):
    """Show the average and confidence interval for HPO terms for a selected population
    Parameters: HPO_list (list): list of HPO # for selected population
                HPO_list_pos (list): list of positive HPO # for selected population
                HPO_list_neg (list): list of negative HPO # for selected population
                logger (Logger): logger
    Returns: None
    Shows the average and CI 95% for HPO counts
    """
    logger.info("HPO pos average : ",np.average(HPO_list_pos),", CI 95% : ",get_CI(HPO_list_pos),", HPO pos max : ",np.max(HPO_list_pos))
    logger.info("HPO neg average : ",np.average(HPO_list_neg),", CI 95% : ",get_CI(HPO_list_neg),", HPO neg max : ",np.max(HPO_list_neg))
    logger.info("HPO average : ",np.average(HPO_list),", CI 95% : ",get_CI(HPO_list),", HPO max : ",np.max(HPO_list))

def get_large_group_HPO(phenotypes):
    """Gets the large groups of HPO hierarchy
    Parameters: phenotypes (pandas.core.DataFrame): pandas dataframe with the phenotypes

    Returns: large_groups_HPO (np.ndarray): list of the large groups in the HPO hierarchy
    """
    large_groups_HPO=[]
    header_phen=list(phenotypes)[1:]
    for phen in header_phen:
        if not(phen.split("\\")[4] in large_groups_HPO):
            large_groups_HPO.append(phen.split("\\")[4])
    return np.array(large_groups_HPO)

def get_phen_to_lg(phenotypes):
    """Gets the large groups of HPO hierarchy
    Parameters: phenotypes (pandas.core.DataFrame): pandas dataframe with the phenotypes

    Returns: list_phenotypes_unique (dict): dict with phenotype as key and large group as value
    """
    header_phen=list(phenotypes)[1:]
    list_phenotypes_unique={}
    for phen in header_phen:
        if not(phen.split("\\")[-2] in list_phenotypes_unique):
            list_phenotypes_unique[phen.split("\\")[-2]]=[phen.split("\\")[4]]
        else:
            if not(phen.split("\\")[4] in list_phenotypes_unique[phen.split("\\")[-2]]):
                list_phenotypes_unique[phen.split("\\")[-2]].append(phen.split("\\")[4])
    return list_phenotypes_unique

def get_large_groups_HPO_count(list_phenotypes_unique, large_groups_HPO, patient_phen, list_patients):
    """Returns the count of HPO terms that belong to a certain group of HPO terms
    Parameters: list_phenotypes_unique (dict): dict with phenotype as key and large group as value
                large_groups_HPO (list): list of large groups that belong to the HPO hierarchy
                patient_phen (dict) : dictionary with patients as keys, with values \
                                    dictionaries with keys ("pos","neg") \
                             with a list of the positive and negative phenotypes presented by each patient
                list_patients (np.ndarray): list of patients to consider
    Returns : group_count : dictionary with keys ("pos","neg") that counts the occurrences of positive or 
                            negative HPO terms for each large group
    """
    group_count={"pos":{lg: 0 for lg in large_groups_HPO},"neg": {lg: 0 for lg in large_groups_HPO}}
    for patient in list_patients:
        for phen in patient_phen[patient]["pos"]:
            for lg in list_phenotypes_unique[phen]:
                group_count["pos"][lg]+=1
        for phen in patient_phen[patient]["neg"]:
            for lg in list_phenotypes_unique[phen]:
                group_count["neg"][lg]+=1
    return group_count
