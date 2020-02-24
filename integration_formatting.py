#!/usr/bin/python
#AUTHOR: RATTINA Vimel - 2020/02/07 - SIB & Enyo Pharma

import sys #I/O files
import os, errno #create folder
import csv #parse tsv
import logging #log file

logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', filename='main.log', level=logging.DEBUG)

###The aim of this script is to create one output file with PPi to be integrated and one output file with interaction mapping information to be integrated

## INPUTS: ##
##-The third-version file (3/4) with PPi and propagated PPi thanks to identical protein sequence

## OUTPUTS: ##
##-An output file with PPi information
##-An output file with interaction mapping information

###################
# FUNCTIONS #######
###################

##Function to store the PPi
def store_ppi(ppi, ppi_list):
    output_line = ppi["pmid"]+"\t"+ppi["psimi_id"]
    output_line += "\t"+ppi["interactor1_accession"]+"\t"+ppi["interactor1_name"]
    output_line += "\t"+ppi["interactor2_accession"]+"\t"+ppi["interactor2_name"]

    if (ppi["ppi_type"] == "vh"):
        output_line += "\t"+ppi["Uniprot_viral_species"]+"\t"+ppi["FTId"]+"\t"+ppi["Chain_name"]
        output_line += "\t"+ppi["interAC_propagation"]
    elif (ppi["ppi_type"] == "hh"):
        output_line += "\t\t\t"
        output_line += "\t"+ppi["interAC_propagation"]

    ppi_list.append(output_line)


##Function to store the interaction mapping information
def store_interaction_mapping(ppi, intmap_list):
    output_line = ppi["pmid"]+"\t"+ppi["psimi_id"]
    output_line += "\t"+ppi["interactor1_accession"]+"\t"+ppi["interactor1_name"]
    output_line += "\t"+ppi["interactor2_accession"]+"\t"+ppi["interactor2_name"]

    if (ppi["interactor1_isoform_accession"] != ""):
        output_line += "\t"+ppi["interactor1_mapping_sequence"]+"\t"+ppi["interactor1_isoform_accession"]
        output_line += "\t"+ppi["interactor1_occurrence_start"]+"\t"+ppi["interactor1_occurrence_stop"]+"\t"+ppi["interactor1_occurrence_identity"]
        
    elif (ppi["interactor2_isoform_accession"] != ""):
        output_line += "\t"+ppi["interactor2_mapping_sequence"]+"\t"+ppi["interactor2_isoform_accession"]
        output_line += "\t"+ppi["interactor2_occurrence_start"]+"\t"+ppi["interactor2_occurrence_stop"]+"\t"+ppi["interactor2_occurrence_identity"]

    if (ppi["ppi_type"] == "vh"):
        output_line += "\t"+ppi["Uniprot_viral_species"]+"\t"+ppi["FTId"]+"\t"+ppi["Chain_name"]
        output_line += "\t"+ppi["interAC_propagation"]
    elif (ppi["ppi_type"] == "hh"):
        output_line += "\t\t\t"
        output_line += "\t"+ppi["interAC_propagation"]

    intmap_list.append(output_line)


##Function to store every interaction mapping information
def integration_formatting(propagated_ppi_file, ppi_output, interaction_mapping_output):
    logging.info("Starts integration formatting")
    ppi_result = open(ppi_output,"w")
    intmap_result = open(interaction_mapping_output, "w")

    ##write the header in PPi output file
    ppi_header = "pmid\tpsimi_id"
    ppi_header += "\tinteractor1_accession\tinteractor1_name"
    ppi_header += "\tinteractor2_accession\tinteractor2_name"
    ppi_header += "\tUniprot_viral_species\tFTId\tChain_name"
    ppi_header += "\tinterAC_propagation"
    ppi_result.write(ppi_header+"\n")
    ##write the header in interaction mapping output file
    intmap_header = "pmid\tpsimi_id"
    intmap_header += "\tinteractor1_accession\tinteractor1_name"
    intmap_header += "\tinteractor2_accession\tinteractor2_name"
    intmap_header += "\tinteractor_mapping_sequence\tinteractor_isoform_accession"
    intmap_header += "\tinteractor_occurrence_start\tinteractor_occurrence_stop\tinteractor_occurrence_identity" 
    intmap_header += "\tUniprot_viral_species\tFTId\tChain_name"
    intmap_header += "\tinterAC_propagation"
    intmap_result.write(intmap_header+"\n")
    
    ppi_list=[]
    intmap_list=[]

    with open(propagated_ppi_file) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        ##fill lists of ppi and interaction mapping
        for ppi in reader:
            store_ppi(ppi, ppi_list)
            if ( (ppi["interactor1_isoform_accession"] != "") or (ppi["interactor2_isoform_accession"] != "") ):
                store_interaction_mapping(ppi, intmap_list)

        ppi_list = list(set(ppi_list)) ##save one copy of the duplicated PPi
        for uniq_ppi in ppi_list:
            ppi_result.write(uniq_ppi+"\n")

        intmap_list = list(set(intmap_list)) ##save one copy of homodimer self-interaction (interactor1 and interactor2 have both the same intmap region)
        for uniq_ppi in intmap_list:
            intmap_result.write(uniq_ppi+"\n")
        
    ppi_result.close()
    intmap_result.close()
    logging.info("Integration formatting ends")

####################
# MAIN #############
####################

if len(sys.argv) < 4:
    sys.exit('Usage: %s <propagated_ppi_file> <ppi_output_file> <interaction_mapping_output_file>\n<propagated_ppi_file>: third version of the file to be integrated in neXtProt\n<ppi_output_file>: PPi file to be integrated in neXtProt\n<interaction_mapping_output_file>: interaction mapping file to be integrated in neXtProt' % sys.argv[0])

if __name__ == "__main__":
    integration_formatting(sys.argv[1], sys.argv[2], sys.argv[3])

#./integration_formatting.py one_iso_integre_vinland-hh-2020-01-22.tsv ppi_integre_vinland-hh-2020-01-22.tsv intmap_integre_vinland-hh-2020-01-22.tsv