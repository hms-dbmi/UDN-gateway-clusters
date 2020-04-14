#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is to conduct the main analysis of the UDN database:
 * Download data using PIC SURE API
 * Analysis of UDN database breaking down patients in adult and pediatric population
 both diagnosed and undiagnosed
 * Clustering of adult and pediatric network using phenotypic similarity with Louvain method
 * Analysis of clusters from a statistical and a phenotypic standpoint
 * Disease enrichment analysis using Orphanet database

example usage from CLI:
 $ python Data_analysis_UDN.py --token personal_token 
                            --json_file "path/to/file" 
                            --genes_file "path/to/gene/info" 
                            -- variants_file "path/to/variant/info" 

For help, run:
 $ Data_analysis_UDN.py -h

"""

__author__ = "Josephine Yates"
__email__ = "josephine.yates@yahoo.fr"

# # Data analysis of UDN patients
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
from scipy.stats import mannwhitneyu, chisquare
from sklearn.metrics.pairwise import pairwise_distances
from docx import Document
from docx.shared import Inches
import ast
import logging

import PicSureHpdsLib
import PicSureClient


# ### Data 

def download_data(resource, logger):
    """download patient information with the PIC SURE API
    Parameters: resource: picsure api resource
                logger (Logger): logger
    Returns: phenotypes, status, genes, variants, primary_symptoms, clinical_site, 
                    family_history, natal_history, demographics, diagnostics (pandas.core.DataFrame):
                    dataframes containing the information for UDN patients
    """

    # download patient data from the server
    phenotypes = get_data_df("\\04_Clinical symptoms and physical findings (in HPO, from PhenoTips)\\", resource)
    status = get_data_df("\\13_Status\\", resource)
    genes=get_data_df("\\11_Candidate genes\\", resource)
    variants=get_data_df("\\12_Candidate variants\\", resource)
    primary_symptoms=get_data_df("\\01_Primary symptom category reported by patient or caregiver\\", resource)
    clinical_site=get_data_df('\\03_UDN Clinical Site\\', resource)
    family_history=get_data_df("\\08_Family history (from PhenoTips)\\", resource)
    natal_history=get_data_df("\\09_Prenatal and perinatal history (from PhenoTips)\\", resource)
    demographics=get_data_df("\\00_Demographics\\", resource)
    diagnostics=get_data_df('\\14_Disorders (in OMIM, from PhenoTips)\\', resource)

    # select only the phenotypes, and not the prenatal phenotypes
    columns_to_del=[]
    for col in list(phenotypes.columns)[1:]:
        if "Prenatal Phenotype" in col.split('\\'):
            columns_to_del.append(col)
    phenotypes=phenotypes.drop(columns_to_del,axis=1)

    # information not in the gateway, example update
    #demographics["\\00_Demographics\\Age at symptom onset in years\\"].loc[patient_not_in_gateway]=0.3
    
    return phenotypes, status, genes, variants, primary_symptoms, clinical_site, \
                    family_history, natal_history, demographics, diagnostics

# ### HPO analysis

def phenotype_formatting(phenotypes, jsonfile, logger):
    # get a dictionnary of patients with their positively and negatively associated phenotypes
    patient_phen = get_patient_phenotypes(phenotypes)

    # retrieve the evaluation date from patients, and remove the information for patients not in the update 
    # (eval date comes from JSON) to fill in 
    patient_eval_date = get_patient_eval_date(jsonfile, patient_phen)

    # get the list of patients evaluated  before 2015, and delete the negative terms for these patients (cf. paper on possible bias)
    # in the entry of negative terms before 2015
    patient_phen = patient_eval_before_2015(patient_eval_date, patient_phen)

    # delete negative terms for patients with over 50 negative terms (cut-off for bias)
    for pat in patient_phen:
        if len(patient_phen[pat]["neg"])>=50:
            patient_phen[pat]["neg"]=[]

    return patient_phen

def patient_breakdown(patient_phen, demographics, status, phenotypes, logger):
    # ### Breakdown into pediatrics, diagnosed or undiagnosed, and adults, diagnosed or undiagnosed
    all_patients = list(patient_phen.keys())
    age = "\\00_Demographics\\Age at symptom onset in years\\"
    adult_patients = \
        demographics[age][demographics[age]>=18.0].index.to_numpy()
    pediatric_patients = \
        demographics[age][demographics[age]<18.0].index.to_numpy()

    # get the list of diagnosed and undiagnosed patients
    list_diagnosed = status.loc[status["\\13_Status\\"] == "solved"].index.values.to_numpy()
    list_undiagnosed = status.loc[status["\\13_Status\\"] != "solved"].index.values.to_numpy()

    # get the list of diagnosed or undiagnosed patients that have at least one HPO term
    list_diagnosed_phen=[patient for patient in list_diagnosed if(patient in patient_phen)]
    list_undiagnosed_phen=[patient for patient in list_undiagnosed if(patient in patient_phen)]


    # get the lists for breakdown (adults diag and undiag, pediatric diag and undiag)
    list_adult_diagnosed=[patient for patient in patient_phen if \
            (patient in adult_patients) and (patient in list_diagnosed_phen)]
    list_adult_undiagnosed=[patient for patient in patient_phen if \
            (patient in adult_patients) and (patient in list_undiagnosed_phen)]
    list_pediatric_diagnosed=[patient for patient in patient_phen if \
            (patient in pediatric_patients) and (patient in list_diagnosed_phen)]
    list_pediatric_undiagnosed=[patient for patient in patient_phen if \
            (patient in pediatric_patients) and (patient in list_undiagnosed_phen)]

    HPO_terms = get_HPO_terms(patient_phen, all_patients)
    HPO_list_pos, HPO_list_neg, HPO_list = get_HPO_count_list(patient_phen, all_patients)
    HPO_list_pos_adult_diagnosed,HPO_list_neg_adult_diagnosed,HPO_list_adult_diagnosed = \
        get_HPO_count_list(patient_phen, list_adult_diagnosed)
    HPO_list_pos_adult_undiagnosed,HPO_list_neg_adult_undiagnosed,HPO_list_adult_undiagnosed = \
        get_HPO_count_list(patient_phen, list_adult_undiagnosed)
    HPO_list_pos_pediatric_diagnosed,HPO_list_neg_pediatric_diagnosed,HPO_list_pediatric_diagnosed = \
        get_HPO_count_list(patient_phen, list_pediatric_diagnosed)
    HPO_list_pos_pediatric_undiagnosed,HPO_list_neg_pediatric_undiagnosed,HPO_list_pediatric_undiagnosed = \
        get_HPO_count_list(patient_phen, list_pediatric_undiagnosed)
    # log the stats on the HPO counts
    show_stats_HPO_counts(HPO_list,HPO_list_pos,HPO_list_neg, logger)
    show_stats_HPO_counts(HPO_list_adult_diagnosed, \
        HPO_list_pos_adult_diagnosed,HPO_list_neg_adult_diagnosed, logger)
    show_stats_HPO_counts(HPO_list_adult_undiagnosed, \
        HPO_list_pos_adult_undiagnosed,HPO_list_neg_adult_undiagnosed, logger)
    show_stats_HPO_counts(HPO_list_pediatric_diagnosed, \
        HPO_list_pos_pediatric_diagnosed,HPO_list_neg_pediatric_diagnosed, logger)
    show_stats_HPO_counts(HPO_list_pediatric_undiagnosed,HPO_list_pos_pediatric_undiagnosed, \
        HPO_list_neg_pediatric_undiagnosed, logger)

    # plot HPO term distribution
    logger.info("Plotting HPO distribution, for all, positive and negative terms.")
    show_distrib_HPO(HPO_list,"Distribution of HPO terms")
    show_distrib_HPO(HPO_list_neg,"Distribution of negative HPO terms")
    show_distrib_HPO(HPO_list_pos,"Distribution of positive HPO terms")

    # get the dataframes with the phenotypes of diagnosed or undiagnosed patients that have at least one HPO term
    phenotypes_diagnosed=phenotypes.loc[list_diagnosed_phen]
    phenotypes_undiagnosed=phenotypes.loc[list_undiagnosed_phen]

    return phenotypes_diagnosed, phenotypes_undiagnosed, list_diagnosed_phen, list_undiagnosed_phen, \
            list_adult_diagnosed, list_adult_undiagnosed, list_pediatric_diagnosed, list_pediatric_undiagnosed, \
            adult_patients, pediatric_patients, HPO_terms, all_patients

# ### HPO large group stats

# get the list of large groups in the HPO hierarchy
def HPO_large_group_analysis(phenotypes, patient_phen, adult_patients, pediatric_patients, all_patients, logger):

    large_groups_HPO = get_large_group_HPO(phenotypes)

    # get the association between unique phenotypes and the large groups \
    # they are related to in the HPO hierarchy
    # list_phenotypes_unique is a dictionary with the phenotypes as keys, \
    # and a list of associated large groups as value
    list_phenotypes_unique = get_phen_to_lg(phenotypes)
    logger.info("Getting count of large groups")
    # get the HPO occurrences for all patients
    large_groups_HPO_count=get_large_groups_HPO_count(list_phenotypes_unique, \
        large_groups_HPO,patient_phen,all_patients)

    logger.info("Total : neg : ",np.sum(list(large_groups_HPO_count["neg"].values())), \
        " pos : ",np.sum(list(large_groups_HPO_count["pos"].values())))

    # get the count of large groups for positive and negative terms of adult patients
    large_groups_HPO_count_adult=get_large_groups_HPO_count(list_phenotypes_unique, \
        large_groups_HPO,patient_phen,adult_patients)

    logger.info("Total adult : neg : ",np.sum(list(large_groups_HPO_count_adult["neg"].values())), \
                " pos : ",np.sum(list(large_groups_HPO_count_adult["pos"].values())))

    # get the count of large groups for positive and negative terms of pediatric patients
    large_groups_HPO_count_pediatric=get_large_groups_HPO_count(list_phenotypes_unique, \
        large_groups_HPO,patient_phen,pediatric_patients)

    logger.info("Total pediatric : neg : ",np.sum(list(large_groups_HPO_count_pediatric["neg"].values())), \
                " pos : ",np.sum(list(large_groups_HPO_count_pediatric["pos"].values())))

    return list_phenotypes_unique, large_groups_HPO

# ### Comparison HPO and Primary Symptoms
def HPO_and_PS(patient_phen, primary_symptoms, list_phenotypes_unique, logger):
    # get the links between the primary symptoms and the HPO large groups
    link_PS_HPO=get_link_between_PS_HPO(patient_phen,primary_symptoms,list_phenotypes_unique)
    return link_PS_HPO

# ### Analysis of metadata
def stats_metadata(demographics, clinical_site, primary_symptoms,
                                    family_history, natal_history,
                                    patient_phen, list_adult_diagnosed, list_adult_undiagnosed, 
                                    list_pediatric_diagnosed, list_pediatric_undiagnosed,
                                    list_diagnosed_phen, list_undiagnosed_phen,
                                    logger):
    # get the dataframes for patients with at least one phenotype, for adult or pediatric, \
    # diagnosed and undiagnosed 
    demographics = demographics.loc[list(patient_phen)]
    demographics_adult_diagnosed = demographics.loc[list_adult_diagnosed]
    demographics_adult_undiagnosed = demographics.loc[list_adult_undiagnosed]
    demographics_pediatric_diagnosed = demographics.loc[list_pediatric_diagnosed]
    demographics_pediatric_undiagnosed = demographics.loc[list_pediatric_undiagnosed]

    clinical_site = clinical_site.loc[list(patient_phen)]
    clinical_site_adult_diagnosed = clinical_site.loc[list_adult_diagnosed]
    clinical_site_adult_undiagnosed = clinical_site.loc[list_adult_undiagnosed]
    clinical_site_pediatric_diagnosed = clinical_site.loc[list_pediatric_diagnosed]
    clinical_site_pediatric_undiagnosed = clinical_site.loc[list_pediatric_undiagnosed]

    cscount_ad = clinical_site_adult_diagnosed.groupby('\\03_UDN Clinical Site\\')['Patient ID'].nunique()
    cscount_and = clinical_site_adult_undiagnosed.groupby('\\03_UDN Clinical Site\\')['Patient ID'].nunique()
    cscount_pd = clinical_site_pediatric_diagnosed.groupby('\\03_UDN Clinical Site\\')['Patient ID'].nunique()
    cscount_pnd = clinical_site_pediatric_undiagnosed.groupby('\\03_UDN Clinical Site\\')['Patient ID'].nunique()

    primary_symptoms = primary_symptoms.loc[list(patient_phen)]
    primary_symptoms_ad = primary_symptoms.loc[list_adult_diagnosed]
    primary_symptoms_and = primary_symptoms.loc[list_adult_undiagnosed]
    primary_symptoms_pd = primary_symptoms.loc[list_pediatric_diagnosed]
    primary_symptoms_pnd = primary_symptoms.loc[list_pediatric_undiagnosed]

    pscount_ad = primary_symptoms_ad.groupby(
        "\\01_Primary symptom category reported by patient or caregiver\\")['Patient ID'].nunique()
    pscount_and = primary_symptoms_and.groupby(
        "\\01_Primary symptom category reported by patient or caregiver\\")['Patient ID'].nunique()
    pscount_pd = primary_symptoms_pd.groupby(
        "\\01_Primary symptom category reported by patient or caregiver\\")['Patient ID'].nunique()
    pscount_pnd = primary_symptoms_pnd.groupby(
        "\\01_Primary symptom category reported by patient or caregiver\\")['Patient ID'].nunique()

    family_history = family_history.loc[list(patient_phen)]
    family_history_ad = family_history.loc[list_adult_diagnosed]
    family_history_and = family_history.loc[list_adult_undiagnosed]
    family_history_pd = family_history.loc[list_pediatric_diagnosed]
    family_history_pnd = family_history.loc[list_pediatric_undiagnosed]

    natal_history = natal_history.loc[list(patient_phen)]
    natal_history = natal_history.replace(0, np.NaN)
    natal_history_ad = natal_history.loc[list_adult_diagnosed]
    natal_history_and = natal_history.loc[list_adult_undiagnosed]
    natal_history_pd = natal_history.loc[list_pediatric_diagnosed]
    natal_history_pnd = natal_history.loc[list_pediatric_undiagnosed]

    # get the statistics for demographics for patient breakdown
    logger.info("Getting the demographics statistics.")
    logger.info("For all patients : {}".format(demographics.describe()))
    logger.info("For adults diagnosed : {}".format(demographics_adult_diagnosed.describe()))
    logger.info("For adults undiagnosed : {}".format(demographics_adult_undiagnosed.describe()))
    logger.info("For pediatrics diagnosed: {}".format(demographics_pediatric_diagnosed.describe()))
    logger.info("For pediatrics undiagnosed : {}".format(demographics_pediatric_undiagnosed.describe()))
    logger.info("Getting the demographics statistics.")
    # get the statistics for clinical sites for patient breakdown
    logger.info("Getting the clinical site statistics.")
    logger.info("For all patients : {}".format(clinical_site.describe()))
    logger.info("For adults diagnosed : {}".format(clinical_site_adult_diagnosed.describe()))
    logger.info("For adults undiagnosed : {}".format(clinical_site_adult_undiagnosed.describe()))
    logger.info("For pediatrics diagnosed: {}".format(clinical_site_pediatric_diagnosed.describe()))
    logger.info("For pediatrics undiagnosed : {}".format(clinical_site_pediatric_undiagnosed.describe()))

    ## can also be done for family and natal history

    natalhist = "\\09_Prenatal and perinatal history (from PhenoTips)\\Maternal Age\\"
    get_diff_parent_age(natal_history_ad, natal_history_pd, natalhist, logger)
    natalhist = "\\09_Prenatal and perinatal history (from PhenoTips)\\Paternal Age\\"
    get_diff_parent_age(natal_history_ad, natal_history_pd, natalhist, logger)

    # show age distribution 
    show_age_distrib(demographics)

    
    # show general statistics for demographics, clinical sites, etc...
    ageeval, ageonset = "\\00_Demographics\\Age at UDN Evaluation (in years)\\", \
                        "\\00_Demographics\\Age at symptom onset in years\\"

    logger.info("Age at UDN Evaluation, adult : {}".format(
        mannwhitneyu(
            np.array(demographics_adult_diagnosed[ageeval]),
            np.array(demographics_adult_undiagnosed[ageeval]))
            ))

    logger.info("Age at UDN Evaluation, pediatric : {}".format(
        mannwhitneyu(
            np.array(demographics_pediatric_diagnosed[ageeval]),
            np.array(demographics_pediatric_undiagnosed[ageeval]))
            ))

    logger.info("Age at symptom onset, adult : {}".format(
        mannwhitneyu(
            np.array(demographics_adult_diagnosed[ageonset]),
            np.array(demographics_adult_undiagnosed[ageonset]))
            ))

    logger.info("Age at symptom onset, pediatric : {}".format(
        mannwhitneyu(
            np.array(demographics_pediatric_diagnosed[ageonset]),
            np.array(demographics_pediatric_undiagnosed[ageonset]))
            ))

    logger.info("Primary symptoms, adults: {}".format(
        mannwhitneyu(np.multiply(list(pscount_ad), 1/len(list_diagnosed_phen)*100),
                    np.multiply(list(pscount_and),1/len(list_undiagnosed_phen)*100))
            ))

    logger.info("Primary symptoms, pediatric {}".format(
        mannwhitneyu(np.multiply(list(pscount_pd), 1/len(list_diagnosed_phen)*100),
                    np.multiply(list(pscount_pnd),1/len(list_undiagnosed_phen)*100))
            ))

    logger.info("Clinical sites, adults: {}".format(
        mannwhitneyu(np.multiply(list(cscount_ad), 1/len(list_diagnosed_phen)*100),
                    np.multiply(list(cscount_and),1/len(list_undiagnosed_phen)*100))
            ))
    
    logger.info("Clinical sites, pediatric {}".format(
        mannwhitneyu(np.multiply(list(cscount_pd), 1/len(list_diagnosed_phen)*100),
                    np.multiply(list(cscount_pnd),1/len(list_undiagnosed_phen)*100))
            ))

    logger.info("Female:male ratio adult diag vs undiag : {}".format(fisher_exact([[22,23],[102,85]])))
    logger.info("Female:male ratio pediatric diag vs undiag : {}".format(fisher_exact([[113,81],[295,321]])))

    return demographics, clinical_site, primary_symptoms, family_history, natal_history

def parental_age_analysis(demographics, natal_history, logger):
    matnatal, patnatal = \
        "\\09_Prenatal and perinatal history (from PhenoTips)\\Maternal Age\\", \
        "\\09_Prenatal and perinatal history (from PhenoTips)\\Paternal Age\\"

    # mat_age is the maternal age without the NaN values
    mat_age=np.array(natal_history[matnatal])
    isnan_mat=np.isnan(mat_age)
    mat_age=mat_age[[not(isnan_mat[i]) for i in range(len(isnan_mat))]]

    # pat_age is the paternal age without the NaN values
    pat_age=np.array(natal_history[patnatal])
    isnan_pat=np.isnan(pat_age)
    pat_age=pat_age[[not(isnan_pat[i]) for i in range(len(isnan_pat))]]

    # distribution of paternal age in the US in 2009 (cf. article)
    USA_dist_pat=[4.7,17.7,25.1,26.6,16.3,6.7,2.1,0.8]
    tranches_pat=["0-19","20-24","25-29","30-34","35-39","40-44","44-50",">50"]
    boundaries_pat=[[0,19],[20,24],[25,29],[30,34],[35,39],[40,44],[44,50],[50,100]]

    # distribution of maternal age in the US in 2009 (cf. article)
    USA_dist_mat = [3.1,6.9,24.4,28.2,23.1,11.5,2.8]
    tranches_mat=["0-18","19","20-24","25-29","30-34","35-39",">39"]
    boundaries_mat=[[0,18],[19,19],[20,24],[25,29],[30,34],[35,39],[39,100]]

    # shows maternal age distribution as opposed to average US in 2009
    dist_age_mat=distrib_age(mat_age,USA_dist_mat,tranches_mat,boundaries_mat,"Maternal")
    ttest_ind(dist_age_mat,USA_dist_mat)

    # shows paternal age distribution as opposed to average US in 2009
    dist_age_pat=distrib_age(pat_age,USA_dist_pat,tranches_pat,boundaries_pat,"Paternal")
    ttest_ind(dist_age_pat,USA_dist_pat)

    # analysis for young versus old paternal age at birth
    paternal_age_df=natal_history[patnatal].dropna()
    patients_with_pat_age=list(paternal_age_df.index)
    old_pat_age=[patients_with_pat_age[i] for i in range(len(patients_with_pat_age)) if \
        paternal_age_df.loc[patients_with_pat_age[i]]>35]
    young_pat_age=[patients_with_pat_age[i] for i in range(len(patients_with_pat_age)) if \
        paternal_age_df.loc[patients_with_pat_age[i]]<=35]

    # Mann Whitney U test for diff between young and old paternal age for Age at UDN evaluation
    ageeval, ageonset = "\\00_Demographics\\Age at UDN Evaluation (in years)\\", \
                        "\\00_Demographics\\Age at symptom onset in years\\"
    mannwhitneyu(list(demographics[ageeval].loc[young_pat_age]),
                    list(demographics[ageeval].loc[old_pat_age]))
    # Mann Whitney U test for diff between young and old paternal age for Age at UDN evaluation
    mannwhitneyu(list(demographics[ageonset].loc[young_pat_age]),
                    list(demographics[ageonset].loc[old_pat_age]))

# ### Genomics
def load_genetic_data(variants_json, genes_json, status, patient_phen, logger):
    variants=get_gene_data(patient_phen, variants_json,"Var")
    genes=get_gene_data(patient_phen, genes_json,"Gene") 

    # get the list of patients that present a candidate gene or candidate variants
    list_patient_genes=list(genes.keys())
    list_patient_variants=list(variants.keys())

    logger.info("Patients in both : {}".format(len([patient for patient in patient_phen if \
        patient in list_patient_genes and patient in list_patient_variants])))
    logger.info("Patients with only genes : {}".format(len([patient for patient in patient_phen if \
        patient in list_patient_genes and not(patient in list_patient_variants)])))
    logger.info("Patients with only variants: {}".format(len([patient for patient in patient_phen if \
        not(patient in list_patient_genes) and patient in list_patient_variants])))

    # count the number of solved cases for people with an indicated gene or an indicated variant
    logger.info("Number of solved and unsolved cases for genes indicated : {}".format( \
                    collec.Counter(status.loc[list(genes.keys())]["\\13_Status\\"])))
    logger.info("Number of solved and unsolved cases for variants indicated : {}".format( \
                    collec.Counter(status.loc[list(variants.keys())]["\\13_Status\\"])))
    return genes, variants

def genomic_analysis(genes, variants, logger):
    # get distribution for variant and gene data
    variant_count=get_dist_genomic(variants,"Var")
    logger.info("Variant distribution : {}".format(variant_count))
    gene_count=get_dist_genomic(genes,"Gen")
    logger.info("Gene distribution : {}".format(gene_count))

    # plot distribution for variants and genes 
    plot_distribution_genomic_data(variants,"Count_dist_var_per_pat.png","variants")
    plot_distribution_genomic_data(genes,"Count_genes_per_pat.png","genes")

# # Clustering

def perform_clustering(phenotypes, patient_phen, adult_patients, pediatric_patients, logger):
    # get the index of unique phenotypes in the phenotype Dataframe
    mat_phen_ind = get_unique_phenotype(phenotypes)
    matrix_phen=phenotypes.drop("Patient ID",axis=1)

    # transform the phenotype dataframe to obtain a matrix of unique phenotypes, \
    # with only patients that have been evaluated,
    # with 1 if the phenotype is positively present, 0 if negative or NaN

    mat_phen_adult=matrix_phen.iloc[:,mat_phen_ind]
    mat_phen_adult=mat_phen_adult.loc[adult_patients]
    mat_phen_adult=mat_phen_adult.replace(to_replace={"Positive": 1, "Negative": 0, np.nan: 0})

    mat_phen_pediatric=matrix_phen.iloc[:,mat_phen_ind]
    mat_phen_pediatric=mat_phen_pediatric.loc[pediatric_patients]
    mat_phen_pediatric=mat_phen_pediatric.replace(to_replace={"Positive": 1, "Negative": 0, np.nan: 0})
    logger.info("Computing jaccard similarity matrix.")
    # we compute the jaccard similarity matrix for the phenotypic matrix, adult patients
    jac_sim_un_adult = 1 - pairwise_distances(mat_phen_adult, metric = "jaccard")

    # we compute the jaccard similarity matrix for the phenotypic matrix, pediatric patients
    jac_sim_un_pediatric = 1 - pairwise_distances(mat_phen_pediatric, metric = "jaccard")

    # create networkx graphs for adult and pediatric network
    logger.info("Creating networkx graphs (long step).")
    # positions can be used to plot the graph in python environment
    graph_un_adult,pos_un_adult=graph_of_patients_js(adult_patients,jac_sim_un_adult)
    graph_un_pediatric,pos_un_pediatric=graph_of_patients_js(pediatric_patients,jac_sim_un_pediatric)

    # writes the computed graph in a gml format, to be able to use Gephi to analyze it further
    logger.info("Writing the graphs in a gml file for analysis with Gephi.")
    nx.write_gml(graph_un_adult,"graph_un_adult.gml")
    nx.write_gml(graph_un_pediatric,"graph_un_pediatric.gml")

    # performing clustering with Louvain method, resolutions 3 and 1.2
    logger.info("Performing clustering (long step).")
    clusters_un_adult=compute_clusters_community(graph_un_adult, 3.0, logger)
    clusters_un_pediatric=compute_clusters_community(graph_un_pediatric, 1.2, logger)

    # get indices of clusters for analysis and outliers
    ind_groups_adult = [cluster for cluster in clusters_un_adult if len(clusters_un_adult[cluster])>3] 
    ind_groups_ped = [cluster for cluster in clusters_un_pediatric if len(clusters_un_pediatric[cluster])>3] 
    ind_outliers_adult = [cluster for cluster in clusters_un_adult if len(clusters_un_adult[cluster])<=3] 
    ind_outliers_ped = [cluster for cluster in clusters_un_pediatric if len(clusters_un_pediatric[cluster])<=3] 

    return clusters_un_adult, clusters_un_pediatric, ind_groups_adult, ind_groups_ped, \
            ind_outliers_adult, ind_outliers_ped

# ### Cluster analysis

def odds_ratio_cluster(clusters_un_adult, clusters_un_pediatric, ind_groups_adult, 
                        ind_groups_pediatric, status, logger):

    OR_diag_adult,IC_adult=calculate_diag_OR(clusters_un_adult,ind_groups_adult,status)
    OR_diag_pediatric,IC_pediatric=calculate_diag_OR(clusters_un_pediatric,ind_groups_pediatric,status)

    return OR_diag_adult, IC_adult, OR_diag_pediatric, IC_pediatric

def phenotype_clusters(clusters_un_adult, clusters_un_pediatric, patient_phen, 
                        ind_groups_adult, ind_groups_pediatric, 
                        ind_outliers_adult, ind_outliers_ped, HPO_terms, logger):

    HPO_count_adult, avg_HPO_clusters_adult, CI_HPO_clusters_adult = \
            get_HPO_count(clusters_un_adult,HPO_terms)
    HPO_count_pediatric, avg_HPO_clusters_pediatric, CI_HPO_clusters_pediatric = \
            get_HPO_count(clusters_un_pediatric,HPO_terms)

    # get the ranked positively and negatively associated phenotypes for patients in each cluster (phen_ranked_pos 
    # and phen_ranked_neg) 
    phen_ranked_pos_adult,phen_ranked_neg_adult = \
            get_phen_ranked(clusters_un_adult, patient_phen, ind_groups_adult)
    phen_ranked_pos_pediatric,phen_ranked_neg_pediatric = \
            get_phen_ranked(clusters_un_pediatric, patient_phen, ind_groups_pediatric)

    # concatenate all outliers for adult and pediatric network
    outlier_pat_adult,outlier_pat_ped = [],[]
    for cl in ind_outliers_adult:
        outlier_pat_adult+=clusters_un_adult[cl]
    for cl in ind_outliers_ped:
        outlier_pat_ped+=clusters_un_pediatric[cl]

    # add the analysis of the outlier population 
    phen_ranked_outliers_adult = phenotype_enrichment_analysis(outlier_pat_adult,patient_phen,"pos")
    phen_ranked_outliers_ped = phenotype_enrichment_analysis(outlier_pat_ped,patient_phen,"pos")
    phen_ranked_pos_adult[len(list(phen_ranked_pos_adult.keys()))]=phen_ranked_outliers_adult
    phen_ranked_pos_pediatric[len(list(phen_ranked_pos_pediatric.keys()))]=phen_ranked_outliers_ped

    # heatmap for positive associations, adult
    logger.info("Phenotype heatmap for adult population")
    heatmap_phen(clusters_un_adult,phen_ranked_pos_adult,ind_groups_adult,
                    "adult",5,12,0,50,"heatmap_adult_clusters")

    # heatmap for positive associations, pediatric
    logger.info("Phenotype heatmap for pediatric population")
    heatmap_phen(clusters_un_pediatric,phen_ranked_pos_pediatric,ind_groups_pediatric,
                    "pediatric",5,12,0,75,"heatmap_pediatric_clusters")

    logger.info("Getting Kruskal Wallis test for HPO count.")
    # get the Kruskal Wallis for the distribution of HPO count between clusters \\ adult
    kr_HPO_ad=kruskal(HPO_count_adult[0],HPO_count_adult[1],HPO_count_adult[2],HPO_count_adult[3])
    logger.info("For adult network : {}".format(kr_HPO_ad))

    # get the Kruskal Wallis for the distribution of HPO count between clusters \\ pediatric
    kr_HPO_ped=kruskal(HPO_count_pediatric[0],HPO_count_pediatric[1],HPO_count_pediatric[2],
                            HPO_count_pediatric[3],HPO_count_pediatric[4])
    logger.info("For pediatric network : {}".format(kr_HPO_ped))

    return phen_ranked_pos_adult, phen_ranked_neg_adult, phen_ranked_pos_pediatric, phen_ranked_neg_pediatric, \
        avg_HPO_clusters_adult, CI_HPO_clusters_adult, avg_HPO_clusters_pediatric, CI_HPO_clusters_pediatric, \
            kr_HPO_ad, kr_HPO_ped, outlier_pat_adult, outlier_pat_ped

def metadata_clusters(clusters_un_adult, clusters_un_pediatric, patient_phen, 
                        ind_groups_adult, ind_groups_pediatric, demographics, logger):

    ageeval, ageonset = '\\00_Demographics\\Age at UDN Evaluation (in years)\\', \
                        '\\00_Demographics\\Age at symptom onset in years\\'
    # get the demographics for the patient in the cluster 
    demographics_coll_adult=metadata_collection(clusters_un_adult,demographics)
    demographics_coll_pediatric=metadata_collection(clusters_un_pediatric,demographics)

    avg_onset_adult,CI_onset_adult = \
        get_metadata_clusters(
            ind_groups_adult,demographics_coll_adult,ageonset,
            )

    avg_UDN_eval_adult,CI_UDN_eval_adult = \
        get_metadata_clusters(
            ind_groups_adult,demographics_coll_adult,ageeval,
            )

    avg_onset_pediatric,CI_onset_pediatric = \
        get_metadata_clusters(
            ind_groups_pediatric,demographics_coll_pediatric,ageonset,
            )

    avg_UDN_eval_pediatric,CI_UDN_eval_pediatric = \
        get_metadata_clusters(
            ind_groups_pediatric,demographics_coll_pediatric,ageeval, 
            )
    # get distributino of gender 
    gender = '\\00_Demographics\\Gender\\'

    gender_distrib_adult = get_distrib(gender,demographics_coll_adult)

    gender_distrib_pediatric = get_distrib(gender,demographics_coll_pediatric)

    logger.info("Getting chi square test to test independence of \
                    female:male ratio depending on cluster presence")
    logger.info(chisquare(np.array([[61,46],[106,102],[21,19],[112,122],[97,106]]).T))

    # KW test for Age at UDN evaluation // adult
    kr_UDN_ad=get_stats_value(ageeval,"adult",ind_groups_adult,demographics_coll_adult, logger)

    # KW test for Age at UDN evaluation // pediatric
    kr_UDN_ped=get_stats_value(ageeval,"pediatric",ind_groups_pediatric,demographics_coll_pediatric, logger)

    # KW test for Age at symptom onset // adult
    kr_onset_ad=get_stats_value(ageonset,"adult",ind_groups_adult,demographics_coll_adult, logger)

    # KW test for Age at symptom onset // pediatric
    kr_onset_ped=get_stats_value(ageonset,"pediatric",ind_groups_pediatric,demographics_coll_pediatric, logger)

    logger.info(
        "Kruskal Wallis test for adult network, age at evaluation : {} and age at symptom onset : {}".format(
            kr_UDN_ad, kr_onset_ad
        )
    )

    logger.info(
        "Kruskal Wallis test for adult network, age at evaluation : {} and age at symptom onset : {}".format(
            kr_UDN_ped, kr_onset_ped
        )
    )
    return avg_onset_adult, CI_onset_adult, avg_onset_pediatric, CI_onset_pediatric, \
            avg_UDN_eval_adult, CI_UDN_eval_adult, avg_UDN_eval_pediatric, CI_UDN_eval_pediatric, \
            kr_UDN_ad, kr_UDN_ped, kr_onset_ad, kr_onset_ped, gender_distrib_adult, gender_distrib_pediatric

# ## Creating tables

def word_tables(clusters_un_adult,avg_HPO_clusters_adult,CI_HPO_clusters_adult,
            gender_distrib_adult,OR_diag_adult,IC_adult,
            avg_onset_adult,CI_onset_adult,avg_UDN_eval_adult,
            CI_UDN_eval_adult, clusters_un_pediatric,avg_HPO_clusters_pediatric,CI_HPO_clusters_pediatric,
            gender_distrib_pediatric,OR_diag_pediatric,IC_pediatric,
            avg_onset_pediatric,CI_onset_pediatric,avg_UDN_eval_pediatric,
            CI_UDN_eval_pediatric, kr_HPO_ad,kr_HPO_ped,kr_UDN_ad,kr_UDN_ped,
            kr_onset_ad,kr_onset_ped, logger):
    
    n_ped,n_ad=5,4


    create_table(n_ad,"A",clusters_un_adult,avg_HPO_clusters_adult,CI_HPO_clusters_adult,
                gender_distrib_adult,OR_diag_adult,IC_adult,
                avg_onset_adult,CI_onset_adult,avg_UDN_eval_adult,
                CI_UDN_eval_adult,"adult_table")

    create_table(n_ped,"P",clusters_un_pediatric,avg_HPO_clusters_pediatric,CI_HPO_clusters_pediatric,
                gender_distrib_pediatric,OR_diag_pediatric,IC_pediatric,
                avg_onset_pediatric,CI_onset_pediatric,avg_UDN_eval_pediatric,
                CI_UDN_eval_pediatric,"pediatric_table")

    create_stat_table(kr_HPO_ad,kr_HPO_ped,kr_UDN_ad,kr_UDN_ped,kr_onset_ad,kr_onset_ped,
                        "statistics adult and pediatric")

# # Diseases enrichment of clusters

def get_mappings(logger):
    # get all the mapping needed from the HPO and Orphadata databases (with additional manual mappings)
    all_diseases = get_all_diseases()
    ORPHAmap, inverseORPHAmap = get_ORPHAmap() 
    mapping_HPO, syn_mapping = get_HPOmap()
    Orphadata_HPO = get_Orphadata_HPO()
    HPO_Orphadata = get_HPO_Orphadata(Orphadata_HPO)
    return all_diseases, ORPHAmap, inverseORPHAmap, mapping_HPO, syn_mapping, Orphadata_HPO, HPO_Orphadata

def get_disease_enrichment_analysis(phen_ranked_pos_adult, phen_ranked_pos_pediatric, 
                                    mapping_HPO, syn_mapping, HPO_Orphadata, 
                                    clusters_un_adult, clusters_un_pediatric,
                                    outlier_pat_adult, outlier_pat_ped, all_diseases, logger):
    
    # get the diseases associated to the most representative phenotypes of the clusters     
    best_phenotypes_dis_analysis_adult = \
        {i: phen_ranked_pos_adult[i][0][:5] for i in list(phen_ranked_pos_adult.keys())}
    best_phenotypes_dis_analysis_ped = \
        {i: phen_ranked_pos_pediatric[i][0][:5] for i in list(phen_ranked_pos_pediatric.keys())}

    # transform into HPO ids
    HPO_ids_cl_adult = get_HPO_from_cluster(best_phenotypes_dis_analysis_adult,mapping_HPO,syn_mapping)
    HPO_ids_cl_ped = get_HPO_from_cluster(best_phenotypes_dis_analysis_ped,mapping_HPO,syn_mapping)

    # get the diseases associated with the HPO ids
    assoc_dis_adult = get_associated_diseases(HPO_ids_cl_adult,HPO_Orphadata)
    assoc_dis_ped = get_associated_diseases(HPO_ids_cl_ped,HPO_Orphadata)

    # concatenate at group level using the Orphadata hierarchy
    assoc_groups_adult = get_associated_groups(assoc_dis_adult,all_diseases)
    assoc_groups_ped = get_associated_groups(assoc_dis_ped,all_diseases)

    # get the relative weight of the disease groups as their size over 
    # the total number of diseases in Orphadata
    tot_weight = np.sum([len(all_diseases[group]) for group in all_diseases])
    group_weight={group: len(all_diseases[group])/tot_weight for group in all_diseases}
    
    # get the weighted # of group assocations 
    weighted_pop_adult = get_weighted_pop(assoc_groups_adult,group_weight)
    weighted_pop_ped = get_weighted_pop(assoc_groups_ped,group_weight)

    # plot as heatmap
    heatmap_real(weighted_pop_adult,clusters_un_adult,len(outlier_pat_adult),
                12,0,15,'heatmap_dis_adult')

    heatmap_real(weighted_pop_ped,clusters_un_pediatric,len(outlier_pat_ped),
                12,0,15,'heatmap_dis_ped')


def main(logger, token):
    # ### Connect to the UDN data resource using the HPDS Adapter
    # Connection to the PicSure Client w/ key
    # token is the individual key given to connect to the resource
    connection = PicSureClient.Client.connect("https://udn.hms.harvard.edu/picsure", token)
    adapter = PicSureHpdsLib.Adapter(connection)
    resource = adapter.useResource("8e8c7ed0-87ea-4342-b8da-f939e46bac26")

    logger.info("Downloading data from server.")
    phenotypes, status, genes, variants, primary_symptoms, clinical_site, \
                    family_history, natal_history, demographics, diagnostics = download_data(resource,logger)
    
    logger.info("Getting the phenotypes for each patient.")
    patient_phen = phenotype_formatting(phenotypes, FLAGS.jsonfile, logger)

    logger.info("Computing analysis for patient breakdown.")
    phenotypes_diagnosed, phenotypes_undiagnosed, list_diagnosed_phen, list_undiagnosed_phen, \
        list_adult_diagnosed, list_adult_undiagnosed, list_pediatric_diagnosed, list_pediatric_undiagnosed, \
        adult_patients, pediatric_patients, HPO_terms, all_patients = \
            patient_breakdown(patient_phen, demographics, status, phenotypes, logger)
    
    logger.info("Computing analysis of HPO large groups.")
    list_phenotypes_unique, large_groups_HPO = HPO_large_group_analysis(phenotypes, 
                        patient_phen, adult_patients, pediatric_patients, all_patients, logger)
    
    logger.info("Getting link between HPO terms and primary symptoms.")
    link_PS_HPO = HPO_and_PS(patient_phen, primary_symptoms, list_phenotypes_unique, logger)

    logger.info("Performing metadata analysis.")
    demographics, clinical_site, primary_symptoms, family_history, natal_history = \
        stats_metadata(demographics, clinical_site, primary_symptoms, family_history, natal_history,
                patient_phen, list_adult_diagnosed, list_adult_undiagnosed, list_pediatric_diagnosed, 
                list_pediatric_undiagnosed, list_diagnosed_phen, list_undiagnosed_phen, logger)
    
    logger.info("Performing parental age analysis.")
    parental_age_analysis(demographics, natal_history, logger)

    logger.info("Downloading genetic information.")
    genes, variants = load_genetic_data(FLAGS.variants_json, FLAGS.genes_json, status, patient_phen, logger)

    logger.info("Performing genomic analysis.")
    genomic_analysis(genes, variants, logger)

    logger.info("Performing clustering on phenotypic data.")
    clusters_un_adult, clusters_un_pediatric, ind_groups_adult, ind_groups_ped, \
            ind_outliers_adult, ind_outliers_ped = perform_clustering(phenotypes, patient_phen, 
                                                                    adult_patients, pediatric_patients, logger)

    logger.info("Starting cluster analysis.")
    # odds ratio 
    OR_diag_adult, IC_adult, OR_diag_pediatric, IC_pediatric = odds_ratio_cluster(clusters_un_adult, 
            clusters_un_pediatric, ind_groups_adult, ind_groups_ped, status, logger)

    # best phenotypes
    phen_ranked_pos_adult, phen_ranked_neg_adult, phen_ranked_pos_pediatric, phen_ranked_neg_pediatric, \
        avg_HPO_clusters_adult, CI_HPO_clusters_adult, avg_HPO_clusters_pediatric, CI_HPO_clusters_pediatric, \
        kr_HPO_ad, kr_HPO_ped, outlier_pat_adult, outlier_pat_ped = \
            phenotype_clusters(clusters_un_adult, clusters_un_pediatric, patient_phen, ind_groups_adult,
                ind_groups_ped, ind_outliers_adult, ind_outliers_ped, HPO_terms, logger)

    # metadata 
    avg_onset_adult, CI_onset_adult, avg_onset_pediatric, CI_onset_pediatric, \
            avg_UDN_eval_adult, CI_UDN_eval_adult, avg_UDN_eval_pediatric, CI_UDN_eval_pediatric, \
            kr_UDN_ad, kr_UDN_ped, kr_onset_ad, kr_onset_ped, \
            gender_distrib_adult, gender_distrib_pediatric = \
            metadata_clusters(clusters_un_adult, clusters_un_pediatric, patient_phen, 
                                ind_groups_adult, ind_groups_ped, demographics, logger)

    logger.info("Saving the analysis tables in word documents.")
    word_tables(clusters_un_adult, avg_HPO_clusters_adult, CI_HPO_clusters_adult, gender_distrib_adult,
                OR_diag_adult, IC_adult, avg_onset_adult, CI_onset_adult, avg_UDN_eval_adult, 
                CI_UDN_eval_adult, clusters_un_pediatric, avg_HPO_clusters_pediatric, 
                CI_HPO_clusters_pediatric, gender_distrib_pediatric, OR_diag_pediatric, 
                IC_pediatric, avg_onset_pediatric, CI_onset_pediatric, avg_UDN_eval_adult,
                CI_UDN_eval_pediatric, kr_HPO_ad, kr_HPO_ped, kr_UDN_ad, kr_UDN_ped,
                kr_onset_ad, kr_onset_ped, logger)

    logger.info("Get mappings from the HPO and the Orphadata databases.") 
    all_diseases, ORPHAmap, inverseORPHAmap, mapping_HPO, syn_mapping, Orphadata_HPO, HPO_Orphadata = \
        get_mappings(logger)
    
    logger.info("Performing diseases enrichment analysis.")
    get_disease_enrichment_analysis(phen_ranked_pos_adult, phen_ranked_pos_pediatric, mapping_HPO,
                                    syn_mapping, HPO_Orphadata, clusters_un_adult, clusters_un_pediatric, 
                                    outlier_pat_adult, outlier_pat_ped, all_diseases, logger)

if __name__=="__main__":

    parser = argparse.ArgumentParser(
        description="CLI args for path to files, individual key for resource"
    )

    parser.add_argument(
        "--token",
        "-t",
        type=str,
        required=True,
        help="personal token granted to access resource",
    )

    parser.add_argument(
        "--json_file",
        "-json",
        type=str,
        required=True,
        help="path to json file containing UDN patient information (might not be necessary if resource updated)",
    )

    parser.add_argument(
        "--genes_file",
        "-genes",
        type=str,
        required=True,
        help="path to json file containing gene information (might not be necessary if resource updated)",
    )

    parser.add_argument(
        "--variants_file",
        "-variants",
        type=str,
        required=True,
        help="path to json file containing variants information (might not be necessary if resource updated)",
    )

    FLAGS = parser.parse_args()

    # clear logger.
    logging.basicConfig(level=logging.DEBUG, filename="log_UDN_analysis.log")

    logger = logging.getLogger("UDN_ANALYSIS")

    # Create a second stream handler for logging to `stderr`, but set
    # its log level to be a little bit smaller such that we only have
    # informative messages
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)

    # Use the default format; since we do not adjust the logger before,
    # this is all right.
    stream_handler.setFormatter(
        logging.Formatter(
            "[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s"
        )
    )
    logger.addHandler(stream_handler)

    logger.info("Usage:\n{0}\n".format(" ".join([x for x in sys.argv])))

    main(logger, FLAGS.token)
