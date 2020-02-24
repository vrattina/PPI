#!/usr/bin/python
#AUTHOR: RATTINA Vimel - 2020/02/11 - SIB & Enyo Pharma

import sys #I/O files
import os, errno #create folder
import csv #parse tsv
import logging #log file
import requests #to POST the sparql request in neXtProt API
import urllib #to urlencode the sparql request

logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', filename='main.log', level=logging.DEBUG)

### This script has DEPENDENCY: /home/vrattina/Documents/HH_PPi/human_identical_proteins.sparql

###The aim of this script is to propagated the PPi if other human proteins were found with similar protein sequence in neXtProt API

## INPUTS: ##
##-The second-version file (2/4) with the interesting information to be integrated in neXtProt and with only one isoform

## OUTPUTS: ##
##-An output file (3/4) with the PPi and propagated PPi thanks to identical protein sequence

###################
# FUNCTIONS #######
###################

##Function sending a request to neXtProt API to retrieve the human proteins with identical amino-acid sequence to other proteins, put everything in a dictionnary
##Dependency: human_identical_proteins.sparql command file and neXtProt API
def find_identical_proteins():
    url = 'https://api.nextprot.org/sparql'
    header = {'accept':'application/sparql-results+json', 'content-type': 'application/x-www-form-urlencoded'}
    query = {'query' : open('/home/vrattina/Documents/HH_PPi/human_identical_proteins.sparql', 'r').read()} #query results in neXtProt are symetric, the interactant will also be the interactor
    payload = urllib.urlencode(query)
    
    r = requests.post(url, payload, headers=header)
    request_dict = r.json()
    identical_prot_dict = {}

    for i in request_dict["results"]["bindings"]:
        identical_proteins = i["identical_protein"]["value"].split(" ")

        protein1 = identical_proteins[0]
        identical_protein_AC1 = protein1.split(":")[0]
        identical_protein_name1 = protein1.split(":")[1]

        protein2 = identical_proteins[1]
        identical_protein_AC2 = protein2.split(":")[0]
        identical_protein_name2 = protein2.split(":")[1]
        
        if identical_prot_dict.has_key(identical_protein_AC1):
            identical_prot_dict[identical_protein_AC1]["partners"].append(identical_protein_AC2)
        else:
            identical_prot_dict[identical_protein_AC1]={}
            identical_prot_dict[identical_protein_AC1]["partners"] = []
            identical_prot_dict[identical_protein_AC1]["partners"].append(identical_protein_AC2)
            identical_prot_dict[identical_protein_AC1]["gene_name"] = identical_protein_name1

    return identical_prot_dict


##Function storing the PPi for the twin protein (having the identical amino-acid sequence) with or without the interaction mapping information
def store_propagated_proteins(identical_prot_dict, current_ppi, interactor_AC, interactor_name, mapping_seq, isoform, occ_start, occ_stop, occ_identity, output_list):
    if ( identical_prot_dict.has_key(interactor_AC) ):
        #if other human protein with identical sequence
        if ( (isoform == "") or (isoform == interactor_AC) ):
            #if the concerned isoform or no interaction mapping information, changes only the AC
            for identical_prot in identical_prot_dict[interactor_AC]["partners"]:
                ppi_twin_protein = current_ppi
                ppi_twin_protein = ppi_twin_protein.replace(interactor_AC, identical_prot)
                ppi_twin_protein = ppi_twin_protein.replace(interactor_name, identical_prot_dict[identical_prot]["gene_name"])
                ppi_twin_protein += "\t"+interactor_AC+"_to_"+identical_prot
                logging.warning(interactor_AC+" get propagated because it has identical amino-acid sequence than "+identical_prot)
                output_list.append(ppi_twin_protein+"\n")
                #print ppi_twin_protein
        # else:
        #     #if not the concerned isoform, save the PPi but not the interaction mapping information
        #     ppi_twin_protein = current_ppi
        #     ppi_twin_protein = ppi_twin_protein.replace(mapping_seq,"").replace(isoform,"")
        #     ppi_twin_protein = ppi_twin_protein.replace(occ_start,"").replace(occ_stop,"")
        #     ppi_twin_protein = ppi_twin_protein.replace(occ_identity,"")

        #     for identical_prot in identical_prot_dict[interactor_AC]["partners"]:
        #         ppi_twin_protein = ppi_twin_protein.replace(interactor_AC, identical_prot)
        #         ppi_twin_protein = ppi_twin_protein.replace(interactor_name, identical_prot_dict[identical_prot]["gene_name"])
        #         ppi_twin_protein += "\t"+interactor_AC+"_to_"+identical_prot+"_VERIFIER"
        #         logging.warning(interactor_AC+" get propagated because it has identical amino-acid sequence than "+identical_prot+" however the interaction mapping is lost because not the corresponding isoform")
        #         #print ppi_twin_protein
        #         output_list.append(ppi_twin_protein+"\n")


##Function writting the output file with the propagated PPi
def inter_ac_propagation(one_isoform_file, propagated_ppi_output):
    logging.info("Starts identical proteins propagation")
    ##retrieve proteins with identical sequence
    identical_prot_dict = find_identical_proteins()
    #print identical_prot_dict

    one_isoform_ppi = open(one_isoform_file, "r")
    propagated_ppi_output = open(propagated_ppi_output, "w")

    header = one_isoform_ppi.readline().rstrip("\n")
    header += "\tinterAC_propagation\n"
    propagated_ppi_output.write(header)

    every_ppi_list = []
    for i in one_isoform_ppi.readlines():
        i = i.rstrip("\n")
        column = i.split("\t")
        ppi_type = column[0]
        interactor_AC1 = column[3]
        interactor_name1 = column[4]
        interactor_AC2 = column[5]
        interactor_name2 = column[6]

        mapping_seq1 = column[7]
        isoform1 = column[8]
        occ_start1 = column[9]
        occ_stop1 = column[10]
        occ_identity1 = column[11]
        
        mapping_seq2 = column[12]
        isoform2 = column[13]
        occ_start2 = column[14]
        occ_stop2 = column[15]
        occ_identity2 = column[16]

        #first store the original ppi
        every_ppi_list.append(i+"\tRaw\n")
        #print i

        #propagate PPi to protein with similar sequence
        if ppi_type == "hh":
            store_propagated_proteins(identical_prot_dict, i,\
                                      interactor_AC1, interactor_name1,\
                                      mapping_seq1, isoform1,\
                                      occ_start1, occ_stop1, occ_identity1,\
                                      every_ppi_list)

            #if non homodimer write it for interactant2
            if ( interactor_AC1 != interactor_AC2 ):
                store_propagated_proteins(identical_prot_dict, i,\
                                          interactor_AC2, interactor_name2,\
                                          mapping_seq2, isoform2,\
                                          occ_start2, occ_stop2, occ_identity2,\
                                          every_ppi_list)
        
        elif ppi_type == "vh":
            #if vh type only interactor_AC1 is the human protein
            #it was decided to do not propagate viral protein
            store_propagated_proteins(identical_prot_dict, i,\
                                    interactor_AC1, interactor_name1,\
                                      mapping_seq1, isoform1,\
                                      occ_start1, occ_stop1, occ_identity1,\
                                      every_ppi_list)
            
    ##some twin proteins are already repertoried by ENYO so no need to put them twice
    #print every_ppi_list
    every_ppi_list = list(set(every_ppi_list))
    for i in every_ppi_list:
        propagated_ppi_output.write(str(i))

    propagated_ppi_output.close()
    one_isoform_ppi.close()
    logging.info("Identical proteins propagation ends")

####################
# MAIN #############
####################

if len(sys.argv) < 3:
    sys.exit('Usage: %s <one_isoform_file> <propagated_ppi_output>\n<one_isoform_file>: second version of the file to be integrated in neXtProt\n<propagated_ppi_output>: third version of the file to be integrated in neXtProt (PPi are propagated to proteins with identical amino-acid sequence)' % sys.argv[0])

if __name__ == "__main__":
    inter_ac_propagation(sys.argv[1], sys.argv[2])

#./inter_ac_propagation.py one_iso_integre_vinland-hh-2020-01-22.tsv toto
