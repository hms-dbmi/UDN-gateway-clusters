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

def get_HPO_from_cluster(best_phenotypes,mapping_HPO,syn_mapping):
    """get the list of HPO ids of the in a cluster using the names 
    Parameters: best_phenotypes (dict): dict with cluster as key, 
                                        list of best phenotypes ranked according to their representation 
                                        in the cluster as value
                mapping_HPO (dict): with HPO terms as value, and dictionnary as value with keys : 
                                    id (the HPO id);
                                    parent (list with parents of the HPO term); 
                                    syn (HPO term synonym if applicable); 
                                    xref (external references numbers, for UMLS, etc.)
                syn_mapping (dict): dictionary that takes into account the synonymy in HPO database, 
                                    with synonymous name as key, and HPO id as value
    Returns: ids (dict): dict with cluster as key, list of associated HPO ids as value
    """
    ids = {cl: [] for cl in best_phenotypes}
    for cl in best_phenotypes:
        for term in best_phenotypes[cl]:
            if term in mapping_HPO:
                ids[cl].append(mapping_HPO[term]["id"])
            else:
                ids[cl].append(syn_mapping[term])
    return ids

def get_associated_diseases(HPO_ids_cl,HPO_Orphadata):
    """get the list of diseases associated with a cluster using the link between HPO and Orphadata 

    Parameters: HPO_ids_cl (dict) : dict with cluster as key, list of associated HPO ids as value
                HPO_Orphadata (dict): dict with Orphanumbers as keys, list of associated HPO as values

    Returns: ass_dis (dict): dict with cluster as key, list of unique associated diseases as values
    """
    ass_dis = {cl: [] for cl in HPO_ids_cl}
    for cl in HPO_ids_cl:
        for HPO in HPO_ids_cl[cl]:
            if HPO in list(HPO_Orphadata.keys()):
                ass_dis[cl]+=HPO_Orphadata[HPO]
        ass_dis[cl]=np.unique(ass_dis[cl])
    return ass_dis

def get_associated_groups(associated_dis,all_diseases):
    """get the count of associations between clusters and the large groups of the HPO hierarchy

    Parameters: associated_dis (dict) : dict with cluster as key, list of unique associated diseases as values
                all_diseases (dict): dict with Orphanet large groups as keys, list of associated diseases as values

    Returns: ass_group (dict): dict with cluster as key, list of associated large groups as values
    """
    ass_group = {cl: [] for cl in associated_dis}
    for cl in associated_dis:
        for dis in associated_dis[cl]:
            for dis_group in list(all_diseases.keys()):
                if dis in all_diseases[dis_group]:
                    ass_group[cl].append(dis_group)
    return ass_group 

def get_weighted_pop(assoc_groups,group_weight):
    """gets the weighted association for each cluster with large groups of Orphadata hierarchy

    Parameters: assoc_group (dict): dict with cluster as key, list of associated large groups as values
                group_weight (dict): dict with group as key, weight as value

    Returns: weighted_pop_cl (dict): dict with cluster as key, dictionary as value with group as key, 
                                group weight as value
    """
    weighted_pop_cl = {cl: {group: 0 for group in group_weight} for cl in assoc_groups}
    for cl in assoc_groups:
        counter = collec.Counter(assoc_groups[cl])
        tot_weighted = np.sum([(1-group_weight[gr])*counter[gr] for gr in list(counter.keys())])
        for gr in list(counter.keys()):
            weighted_pop_cl[cl][gr]=np.log(1+((1-group_weight[gr])*counter[gr])/tot_weighted)*100
    return weighted_pop_cl

### Mappings

def get_list_diseases_from_XML(filename):
    """get the list of diseases in the Orphanet database by disease groups
    Parameters: filename (str): file name for the Orphanet xml file
    Returns: dis (list): list of diseases for the file name 
    """
    tree = ET.parse(filename)
    root = tree.getroot()
    dis=[]
    for s in root.findall(".//Disorder"):
        dis.append(s.find("OrphaNumber").text.lower())
    return dis

def get_all_diseases():
    """get the list of all diseaes associated to their group in the orphanet hierarchy 
    Parameters: None
    Returns: all_diseases (dict): dictionary with hpo groups as keys, list of diseases as value
    """
    all_diseases={}
    for file in os.listdir("xml_orphanet/"):
        if file.endswith(".txt"):
            diseases=get_list_diseases_from_XML("xml_orphanet/"+str(file))
            all_diseases[file.split(".")[0]]=diseases
    return all_diseases

def get_ORPHAmap():
    """get the mapping (and inverse mapping) between orphanet names and orphanumbers
    Parameters: None
    Returns: ORPHAmap (dict): dictionary with orphanumbers as keys, orphanet names as values 
            inverseORPHAmap (dict): dictionary with orphanet names as keys, orphanumbers as values 
    """
    ORPHAmap,inverseORPHAmap={},{}
    for file in os.listdir("xml_orphanet/"):
        if file.endswith(".xml"):
            tree = ET.parse('xml_orphanet/'+file)
            root = tree.getroot()
            for s in root.findall(".//Disorder"):
                orphan=s.find("OrphaNumber").text
                name=s.find("Name").text
                ORPHAmap[orphan]=name
                inverseORPHAmap[name]=orphan
    # manual mapping for missing diseases (incoherence between databases)
    ORPHAmap["2168"]="Homocarnosinosis"
    inverseORPHAmap["Homocarnosinosis"]="2168"
    ORPHAmap["100990"]="Autosomal dominant spastic paraplegia type 9"
    inverseORPHAmap["Autosomal dominant spastic paraplegia type 9"]="100990"
    ORPHAmap["98028"]="Rare circulatory system disease"
    inverseORPHAmap["Rare circulatory system disease"]="98028"
    ORPHAmap["77"]="Aniridia"
    inverseORPHAmap["Aniridia"]="77"
    ORPHAmap["2580"]="Shoulder and girdle defects-familial intellectual disability syndrome"
    inverseORPHAmap["Shoulder and girdle defects-familial intellectual disability syndrome"]="2580"
    ORPHAmap["956"]="Acropectororenal dysplasia"
    inverseORPHAmap["Acropectororenal dysplasia"]="956"
    ORPHAmap["2286"]="Solitary median maxillary central incisor syndrome"
    inverseORPHAmap["Solitary median maxillary central incisor syndrome"]="2286"
    ORPHAmap["3105"]="Robinow-like syndrome"
    inverseORPHAmap["Robinow-like syndrome"]="3105"
    ORPHAmap["1396"]="Cerebrorenodigital syndrome"
    inverseORPHAmap["Cerebrorenodigital syndrome"]="1396"
    ORPHAmap["2894"]="Pilotto syndrome"
    inverseORPHAmap["Pilotto syndrome"]="2894"
    ORPHAmap["1092"]="Renal-genital-middle ear anomalies"
    inverseORPHAmap["Renal-genital-middle ear anomalies"]="1092"
    ORPHAmap["171836"]="Amelogenesis imperfecta-gingival hyperplasia syndrome"
    inverseORPHAmap["Amelogenesis imperfecta-gingival hyperplasia syndrome"]="171836"

    return ORPHAmap, inverseORPHAmap

def get_HPOmap():
    """get the HPO mapping and inverse mapping using the hpo txt file
    Parameters: None
    Returns: mapping_HPO (dict): with HPO terms as value, and dictionnary as value with keys : 
                                    id (the HPO id);
                                    parent (list with parents of the HPO term); 
                                    syn (HPO term synonym if applicable); 
                                    xref (external references numbers, for UMLS, etc.)
            syn_mapping (dict): dictionary that takes into account the synonymy in HPO database, 
                                    with synonymous name as key, and HPO id as value
    """
 
    mapping_HPO={}
    with open("hpo.txt","r+") as hpo:
        lines = hpo.readlines()
        i=0
        for i in range(len(lines)):
            if lines[i].split(":")[0]=="id":
                hpoid=lines[i].split(" ")[1].split("\n")[0]
                name=""
                for namestr in lines[i+1].split(" ")[1:]:
                    name+=namestr+" "
                name=name.split("\n")[0]
                mapping_HPO[name]={"id": hpoid, "xref": [], "parent": [], "syn":[]}
            if lines[i].split(" ")[0]=="xref:":
                mapping_HPO[name]["xref"].append(lines[i].split(" ")[1].split("\n")[0])
            if lines[i].split(" ")[0]=="is_a:":
                mapping_HPO[name]["parent"].append(lines[i].split(" ")[1])
            if lines[i].split(" ")[0]=="synonym:":
                namesyn=""
                for namestr in lines[i].split(" ")[1:]:
                    namesyn+=namestr+" "
                namesyn=namesyn.split("\"")[1].split("'")[0]
                if len(namesyn.split("obsolete "))>1:
                    namesyn=namesyn.split("obsolete ")[1]
                mapping_HPO[name]["syn"].append(namesyn)

    # additionnal manual mapping (for mistakes/missing tokens/missing terms...)
    syn_mapping={}
    for dis in mapping_HPO:
        for syn in mapping_HPO[dis]["syn"]:
            syn_mapping[syn]=mapping_HPO[dis]["id"]
    syn_mapping["Contracture of the distal interphalangeal joints of the fingers"]="HP:0009697"
    syn_mapping["Decreased testosterone in males"]="HP:0040171"
    syn_mapping["EMG myopathic abnormalities"]="HP:0003458"
    syn_mapping["EMG myotonic discharges"]="HP:0100284"
    syn_mapping["Increased IgE level"]="HP:0003212"
    syn_mapping["EMG chronic denervation signs"]="HP:0003444"
    syn_mapping["Primitive reflexes (palmomental, snout, glabellar)"]="HP:0002476"
    syn_mapping["EMG slow motor conduction"]="HP:0100287"
    syn_mapping["Severe Myopia"]="HP:0011003"
    syn_mapping["Arthralgiaarthritis"]="HP:0005059"
    syn_mapping["Decreased CSF homovanillic acid (HVA)"]="HP:0003785"
    syn_mapping["Hip Subluxation"]="HP:0030043"
    syn_mapping["Increased IgM level"]="HP:0003496"
    syn_mapping["Hyperpigmentedhypopigmented macules"]="HP:0007441"
    syn_mapping["Noninflammatory macular atrophy"]="HP:0007401"
    syn_mapping["Capillary hemangiomas"]="HP:0005306"
    syn_mapping["AplasiaHypoplasia of the lungs"]="HP:0006703"
    syn_mapping["Abnormal serum cobalamin"]="HP:0040126"
    syn_mapping["Cone-rod dystrophy"]="HP:0000548"
    syn_mapping["AplasiaHypoplasia of the inner ear"]="HP:0008774"
    syn_mapping["Decreased number of CD4+ T cells"]="HP:0032183"
    syn_mapping["Enlarged kidneys"]="HP:0000113"
    syn_mapping["Nephroblastoma (Wilms tumor)"]="HP:0002667"
    syn_mapping["Cervical vertebral fusion (C2C3)"]="HP:0002949"
    syn_mapping["Pulmonary hypertension"]="HP:0030950"
    syn_mapping["Prominent epicanthal folds"]="HP:0000286"
    syn_mapping["Cellulitis due to immunodeficiency"]="HP:0100658"
    syn_mapping["Abnormal brain cholinecreatine ratio by MRS"]="HP:0012709"
    syn_mapping["Small palpebral fissure"]="HP:0045025"
    syn_mapping["EMG impaired neuromuscular transmission"]="HP:0100285"
    syn_mapping["Absenthypoplastic coccyx"]="HP:0008436"
    syn_mapping["Short tubular bones (hand)"]="HP:0001248"
    syn_mapping["Increased serum Insulin-like growth factor 1"]="HP:0030269"
    syn_mapping["AplasiaHypoplasia of the middle phalanx of the 4th toe"]="HP:0100373"
    syn_mapping["AplasiaHypoplasia of the middle phalanx of the 5th toe"]="HP:0100374"
    syn_mapping["Cortical thickening (humeral)"]="HP:0003868"
    syn_mapping["Abnormality of Descemet's membrane"]="HP:0011490"
    syn_mapping["Increased serum free triiodothyronine (fT3)"]="HP:0011788"
    syn_mapping["Limited knee flexionextension"]="HP:0005085"
    syn_mapping["Limited pronationsupination of forearm"]="HP:0006394"
    syn_mapping["Biilateral vocal cord paralysis"]="HP:0012820"
    syn_mapping["Macroreticular retinal dystrophy"]="HP:0000556"
    syn_mapping["Congenital visual impairment"]="HP:0000505"
    syn_mapping["AtrophyDegeneration affecting the brainstem"]="HP:0007366"
    syn_mapping["Midface prominence"]="HP:0430026"
    syn_mapping["Basal lamina 'onion bulb' formation"]="HP:0003400"
    syn_mapping["Neutrophillia"]="HP:0011897"
    syn_mapping["Elevated circulating parathyroid hormone (PTH) level"]="HP:0003165"
    syn_mapping["Generalized cerebral atrophyhypoplasia"]="HP:0007058"
    syn_mapping["Somnolence"]="HP:0001262"
    syn_mapping["Increased IgA level"]="HP:0003261"
    syn_mapping["Increased blood urea nitrogen (BUN)"]="HP:0003138"
    syn_mapping["Decreased number of CD8+ T cells"]="HP:0410385"
    syn_mapping["Abnormality of acetylcarnitine metabolism"]="HP:0012071"
    syn_mapping["Hemisacrum (S2-S5)"]="HP:0009790"
    syn_mapping["EMG Positive sharp waves"]="HP:0030007"
    syn_mapping["Spastichyperactive bladder"]="HP:0005340"
    syn_mapping["AplasiaHypoplasia of the clavicles"]="HP:0006710"
    syn_mapping["Congenital exotropia"]="HP:0000577"
    syn_mapping["Abnormal rapid eye movement (REM) sleep"]="HP:0002494"
    syn_mapping["Flared metaphyses (elbow)"]="HP:0003950"
    syn_mapping["Cortical subperiosteal resorption (humeral metaphyses)"]="HP:0003909"
    syn_mapping["Abnormality of ornithine metabolism"]="HP:0012025"
    syn_mapping["Oligodactyly (feet)"]="HP:0001849"
    syn_mapping["AplasiaHypoplasia of the musculature of the pelvis"]="HP:0001471"
    syn_mapping["Status Asthmaticus"]="HP:0012653"
    syn_mapping["Vitreoretinal abnormalities"]="HP:0007773"
    syn_mapping["Abnormality of natural killer cell number"]="HP:0040089"
    syn_mapping["AplasiaHypoplasia of the colon"]="HP:0100811"
    syn_mapping["AplasiaHypoplasia of the nasal bone"]="HP:0010940"
    syn_mapping["Abnormality of B cell number"]="HP:0010975"
    syn_mapping["Abnormality of lymphocytes"]="HP:0004332"
    syn_mapping["AplasiaHypoplasia of the eyebrow"]="HP:0100840"
    syn_mapping["AplasiaHypoplasia of the mandible"]="HP:0009118"
    syn_mapping["AplasiaHypoplasia of the distal phalanges of the hand"]="HP:0009835"
    syn_mapping["AplasiaHypoplasia of the middle phalanges of the hand"]="HP:0009843"
    syn_mapping["AplasiaHypoplasia of the thumb"]="HP:0009601"
    syn_mapping["Oligodactyly (hands)"]="HP:0012165"
    syn_mapping["AplasiaHypoplasia of the middle phalanx of the 3rd toe"]="HP:0100372"
    syn_mapping["AplasiaHypoplasia of the 3rd toe"]="HP:0010331"
    syn_mapping["AplasiaHypoplasia of the 4th toe"]="HP:0010337"
    syn_mapping["AplasiaHypoplasia of the 5th toe"]="HP:0010343"
    syn_mapping["AplasiaHypoplasia of toe"]="HP:0001991"
    syn_mapping["Abnormality of the metaphyses"]="HP:0003907"
    syn_mapping["AplasiaHypoplasia of the ulna"]="HP:0006495"
    syn_mapping["Abnormality of the fibula"]="HP:0010595"
    syn_mapping["Abnormality of the metatarsal bones"]="HP:0001832"
    syn_mapping["Abnormality of the tibia"]="HP:0002992"
    syn_mapping["Decreasedabsent ankle reflexes"]="HP:0200101"
    syn_mapping["Abnormal enzymecoenzyme activity"]="HP:0012379"
    syn_mapping["Abnormality of carbohydrate metabolismhomeostasis"]="HP:0011013"
    syn_mapping["Abnormality of aspartate family amino acid metabolism"]="HP:0010899"
    syn_mapping["Abnormality of citrulline metabolism"]="HP:0011965"
    syn_mapping["Abnormality of homocysteine metabolism"]="HP:0010919"
    syn_mapping["Chromsome breakage"]="HP:0040012"
    syn_mapping["Decreased activity of the pyruvate dehydrogenase (PDH) complex"]="HP:0002928"
    syn_mapping["Abnormality of copper homeostasis"]="HP:0010836"
    syn_mapping["Abnormal serum iron"]="HP:0040130"
    syn_mapping["Abnormality of lipid metabolism"]="HP:0003119"
    syn_mapping["Abnormality of fatty-acid metabolism"]="HP:0004359"
    syn_mapping["Abnormality of carnitine metabolism"]="HP:0010967"
    syn_mapping["Abnormality of long-chain fatty-acid metabolism"]="HP:0010964"
    syn_mapping["Abnormality of glycine metabolism"]="HP:0010895"
    syn_mapping["Abnormality of the esophagus"]="HP:0025270"
    syn_mapping["Elevated hepatic transaminases"]="HP:0002910"
    syn_mapping["Abnormality of cardiac atrium"]="HP:0005120"
    syn_mapping["Ebstein's anomaly of the tricuspid valve"]="HP:0010316"
    syn_mapping["Effort-induced polymorphic ventricular tachycardias"]="HP:0004758"
    syn_mapping["Abnormality of circle of Willis"]="HP:0012518"
    syn_mapping["Coronary artery disease"]="HP:0006704"
    syn_mapping["Peripheral arterial disease"]="HP:0004950"
    syn_mapping["Dilatation of the ascending aorta"]=""
    syn_mapping["AplasiaHypoplasia of the tragus"]="HP:0009913"
    syn_mapping["Adrenocorticotropin (ACTH) deficient adrenal insufficiency"]="HP:0011735"
    syn_mapping["Abnormality of the conjunctiva"]="HP:0008054"
    syn_mapping["Congenital primary aphakia"]=""
    syn_mapping["AplasiaHypoplasia of the optic nerve"]="HP:0008058"
    syn_mapping["Abnormality of the vitreous humor"]="HP:0004327"
    syn_mapping["Congenital strabismus"]="HP:0000486"
    syn_mapping["Abnormality of vision evoked potentials"]="HP:0000649"
    syn_mapping["Hemianopic blurring of vision"]="HP:0001125"
    syn_mapping["Congenital glaucoma"]="HP:0008007"
    syn_mapping["AplasiaHypoplasia of the testes"]="HP:0010468"
    syn_mapping["Aplastichypoplastic toenail"]="HP:0010624"
    syn_mapping["AplasiaHypoplasia of the nails"]="HP:0008386"
    syn_mapping["EMG axonal abnormality"]="HP:0003482"
    syn_mapping["EMG neuropathic changes"]="HP:0003445"
    syn_mapping["Abnormality of the globus pallidus"]="HP:0002454"
    syn_mapping["AplasiaHypoplasia of the corpus callosum"]="HP:0007370"
    syn_mapping["AplasiaHypoplasia of the cerebral white matter"]="HP:0012429"
    syn_mapping["AtrophyDegeneration affecting the cerebrum"]="HP:0007369"
    syn_mapping["Porencephaly"]=""
    syn_mapping["AplasiaHypoplasia of the cerebellar vermis"]="HP:0006817"
    syn_mapping["Brain very small"]=""
    syn_mapping["AtrophyDegeneration involving the spinal cord"]="HP:0007344"
    syn_mapping["AplasiaHypoplasia involving the central nervous system"]="HP:0002977"
    syn_mapping["AtrophyDegeneration affecting the central nervous system"]="HP:0007367"
    syn_mapping["CNS infection"]="HP:0011450"
    syn_mapping["Abnormal pyramidal signs"]="HP:0007256"
    syn_mapping["Hemiplegiahemiparesis"]="HP:0004374"
    syn_mapping["Reduced consciousnessconfusion"]="HP:0004372"
    syn_mapping["Inability to walk by childhoodadolescence"]="HP:0006915"
    syn_mapping["Abnormal emotionaffect behavior"]="HP:0100851"
    syn_mapping["Abnormal fearanxiety-related behavior"]="HP:0100852"
    syn_mapping["Abnormality of the lung"]="HP:0002088"
    syn_mapping["Abnormality of the tracheobronchial system"]=""
    syn_mapping["AplasiaHypoplasia of the ribs"]="HP:0006712"
    syn_mapping["Chronic recurrent multifocal osteomyelitis"]=""
    syn_mapping["Increased number of peripheral CD3+ T cells"]=""

    return mapping_HPO, syn_mapping

def get_Orphadata_HPO():
    """get the mapping between Orphadata and HPO values using preprocessed dictionary
    Parameters: None
    Returns: Orphadata_HPO (dict): dict with Orphanumbers as keys, list of associated HPO as values
    """
    with open('xml_orphanet/orphadata_dict.csv') as csv_file:
        reader = csv.reader(csv_file)
        Orphadata_HPO = dict(reader)
    for dis in Orphadata_HPO:
        Orphadata_HPO[dis] = np.array(ast.literal_eval(Orphadata_HPO[dis]))
    return Orphadata_HPO

def get_HPO_Orphadata(Orphadata_HPO):
    """get the mapping between Orphadata and HPO values using preprocessed dictionary
    Parameters: None
    Returns: Orphadata_HPO (dict): dict with Orphanumbers as keys, list of associated HPO as values
    """
    HPO_Orphadata={}
    for dis in list(Orphadata_HPO.keys()):
        for HPO in Orphadata_HPO[dis]:
            if HPO in list(HPO_Orphadata.keys()):
                HPO_Orphadata[HPO].append(dis)
            else:
                HPO_Orphadata[HPO] = [dis]
    return HPO_Orphadata