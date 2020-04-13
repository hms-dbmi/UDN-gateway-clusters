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

def calculate_diag_OR(clusters,clusters_ind,status):
    """Calculate the Odds Ratio for the probability of being diagnosed linked to being in a certain cluster
    Parameters: clusters (dict): dictionary with cluster number as key and list of patients in cluster as value
                clusters_ind (list): indices of cluster to take into account 
                status (str): status of the patient (if patient's case is solved or not)
    Returns: OR_diag (dict): dictionary with cluster number as key and the Odds Ratio (OR) for each cluster
    """
    count_diag_clusters={cluster: 0 for cluster in clusters_ind}
    for cluster in clusters_ind:
        for patient in clusters[cluster]:
            if status.loc[patient]["\\13_Status\\"]=="solved":
                count_diag_clusters[cluster]+=1
    OR_diag,IC={},{}
    def IC_func(sign,OR,a,b,c,d):
        if (a==0 or b==0 or c==0 or d==0):
            return None
        if sign=="up":
            return np.exp(np.log(OR)+1.96*np.sqrt(1/a+1/b+1/c+1/d))
        else:
            return np.exp(np.log(OR)-1.96*np.sqrt(1/a+1/b+1/c+1/d))
    for cluster in count_diag_clusters:
        count_diag_notin=np.sum([count_diag_clusters[cl] for cl in clusters_ind if not(cl==cluster)])
        OR_diag[cluster] = \
            (
                (count_diag_clusters[cluster]/count_diag_notin)/ \
                ((len(clusters[cluster])-count_diag_clusters[cluster])/ \
                                np.sum([len(clusters[cl])-count_diag_clusters[cl] for cl in clusters_ind]))
                    )
        IC[cluster]={"up": IC_func(
                            "up",OR_diag[cluster],count_diag_clusters[cluster],
                            (len(clusters[cluster])-count_diag_clusters[cluster]),
                            count_diag_notin,
                            np.sum([len(clusters[cl])-count_diag_clusters[cl] for cl in clusters_ind]),
                            ),
                    "low": IC_func(
                            "low",OR_diag[cluster],count_diag_clusters[cluster],
                            (len(clusters[cluster])-count_diag_clusters[cluster]),
                            count_diag_notin,
                            np.sum([len(clusters[cl])-count_diag_clusters[cl] for cl in clusters_ind])
                            )
                    }
    return OR_diag,IC

def phenotype_enrichment_analysis(patients_clustered,patient_phen,polarity_HPO):
    """Get the phenotypes shared by the most patients in the cluster according to polarity (positive or negative)
    Parameters: patients_clustered (list): list of patients in the cluster 
                patient_phen (dict): dictionary of unique phenotypes associated with each patient; 
                        key is patient, value is dictionary with key "pos" or "neg" 
                        and value list of unique phenotypes with positive or negative association
                polarity_HPO (str): "pos" or "neg", polarity wanted for the phenotype enrichment analysis

    Returns: phen_ranked (np.ndarray): list of best phenotypes ranked according to their representation 
                                in the cluster
             values (np.ndarray): list of proportion of patients presenting the phenotype in 
                                the phen_ranked same position (ex: values[i] will have the 
                                represention of phenotype phen_ranked[i])
    """
    phen_count={}
    for patient in patients_clustered:
        for phen in patient_phen[patient][polarity_HPO]:
            if not(phen in phen_count):
                phen_count[phen]=1/len(patients_clustered)
            else:
                phen_count[phen]+=1/len(patients_clustered)
    phen_ranked=np.array([phen for phen in phen_count.keys()])
    values=np.array([phen_count[phen] for phen in phen_ranked])
    indrank=np.argsort(values)[::-1]
    phen_ranked=phen_ranked[indrank]
    values=values[indrank]
    return phen_ranked,values

def get_HPO_count(patients_clustered,HPO_terms):
    """get the count of HPO terms for patients in the cluster, and the average
    Parameters: patients_clustered (dict): dictionary with cluster number as key and list 
                                of patients in the cluster as value
                HPO_terms (dict): dictionary with patient as key and 
                            count of HPO terms for the patient as value
    Returns: HPO_cluster (dict): dictionary with cluster number as key and list of HPO numbers 
                            for each patient in the cluster as value
             avg_HPO_clusters (dict): dictionary with cluster number as key and 
                            average number of HPO terms per patient as value
    """
    HPO_cluster = {i: [] for i in patients_clustered.keys()}
    for cluster in patients_clustered:
        for patient in patients_clustered[cluster]:
            HPO_cluster[cluster].append(HPO_terms[patient])
    avg_HPO_clusters = {cluster: np.average(HPO_cluster[cluster]) for cluster in patients_clustered.keys()}
    CI_HPO_clusters = {cluster: get_CI(HPO_cluster[cluster]) for cluster in patients_clustered.keys()}
    return HPO_cluster,avg_HPO_clusters,CI_HPO_clusters

def get_phen_ranked(clusters_un,ind_groups):
    """get the count of HPO terms for patients in the cluster, and the average
    Parameters: clusters_un (dict): dictionary with cluster number as key and list of patients in cluster as value
                ind_groups (list): list of clusters to analyze
    Returns: phen_ranked_pos (dict): dictionary with cluster as key, array of sorted positive 
                            phenotypes and proportion values as value
            phen_ranked_neg (dict): dictionary with cluster as key, array of sorted negative 
                            phenotypes and proportion values as value
    """
    phen_ranked_pos,phen_ranked_neg={cluster: [] for cluster in ind_groups},{cluster: [] for cluster in ind_groups}
    for cluster in ind_groups:
        phen_ranked_pos[cluster]=phenotype_enrichment_analysis(clusters_un[cluster],patient_phen,"pos")
        phen_ranked_neg[cluster]=phenotype_enrichment_analysis(clusters_un[cluster],patient_phen,"neg")
    return phen_ranked_pos,phen_ranked_neg

def metadata_collection(patients_clustered,metadata):
    """Get the metadata for each cluster 
    Parameters: patients_clustered (dict): dictionary with cluster number as key and 
                                    list of patients in the cluster as value
                metadata (pandas.core.DataFrame): dataframe with metadata
    Returns: metadata_clusters (dict): dictionary with clusters as keys and dictionary as value, 
                                    with key the metadata considered and list of values for 
                                    patients in the cluster as value
    """
    metadata_clusters={cl: {meta: [] for meta in list(metadata.columns)} \
                                for cl in patients_clustered.keys()}
    for cl in patients_clustered:
        for patient in patients_clustered[cl]:
            for meta in list(metadata.columns)[1:]:
                metadata_clusters[cl][meta].append(metadata.loc[patient][meta])
    return metadata_clusters

def get_metadata_clusters(ind_groups,demographics_coll,attribute):
    """Returns the average and 95% CI for an attribute
    Parameters: ind_groups (list): indices of clusters to consider
                demographics_coll (list): dictionary with clusters as keys and dictionary as value, 
                                    with key the metadata considered and list of values for patients 
                                    in the cluster as value
                attribute (str): attribute to consider (must be in demographics_coll)

    Returns: avg_att (dict): dictionnary with clusters as keys and average of considered attribute as value
             CI_att (dict): dictionnary with clusters as keys and tuple with lower and upper CI 95% as value"""
    avg_att,CI_att={cl: 0 for cl in ind_groups},{cl: (0,0) for cl in ind_groups}
    for cluster in ind_groups:
        lst=np.array(demographics_coll[cluster][attribute])
        lst=lst[np.logical_not(np.isnan(lst))]
        avg_att[cluster]=np.average(lst)
        CI_att[cluster]=get_CI(lst)
    return avg_att,CI_att

def get_distrib(attribute,demographics_coll):
    """Get a distribution for an attribute
    Parameters: attribute (str): string, attribute we want the distribution of
    Returns: counter (dict): counter of the collection library (distribution)
    """
    counter={}
    for cluster in demographics_coll:
        dc=np.array(demographics_coll[cluster][attribute])
        if type(dc[0])==np.float64:
            dc=dc[np.logical_not(np.isnan(dc))]
        counter[cluster]=dict(collec.Counter(dc))
    return counter

def get_stats_value(value_considered,ad_or_ped,ind_clusters,demographics_coll,logger):
    """get the Kruskal Wallis U index and p-value for a type of demographics
    Parameters: value_considered (str): type of demographics we want the KW U index and p-value
                ad_or_ped (str): "adult" or "pediatric", to consider one or the other cluster
                ind_clusters (list): indices of clusters to consider
                demographics_coll (dict): dictionary with clusters as keys and dictionary as value, 
                                        with key the metadata considered and list of values for 
                                        patients in the cluster as value
                logger (Logger): logger

    Returns: (tuple) Kruskal-Wallis statistic, 2-D tuple with H-index (0) and p-value (1)
    Prints the KW H index and p-value
    """
    dc={i: [] for i in ind_clusters}
    for i in ind_clusters:
        dc[i]=np.array(demographics_coll[i][value_considered])
        if type(dc[i][0])==np.float64:
            dc[i]=dc[i][np.logical_not(np.isnan(dc[i]))]
    if ad_or_ped=="adult":
        logger.info("Kruskal Wallis statistics for adult network.")
        logger.info(kruskal(dc[0],dc[1],dc[2],dc[3]))
        return kruskal(dc[0],dc[1],dc[2],dc[3])
    elif ad_or_ped=="pediatric":
        logger.info("Kruskal Wallis statistics for pediatric network.")
        logger.info(kruskal(dc[0],dc[1],dc[2],dc[3],dc[4]))
        return kruskal(dc[0],dc[1],dc[2],dc[3],dc[4])
