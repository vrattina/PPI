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
    output_line += "\tNX_"+ppi["interactor1_accession"]+"\t"+ppi["interactor1_name"]
    output_line += "\tNX_"+ppi["interactor2_accession"]+"\t"+ppi["interactor2_name"]

    # if (ppi["ppi_type"] == "vh"):
    #     output_line += "\t"+ppi["Uniprot_viral_species"]+"\t"+ppi["FTId"]+"\t"+ppi["Chain_name"]
    #     output_line += "\t"+ppi["interAC_propagation"]
    # elif (ppi["ppi_type"] == "hh"):
    #     output_line += "\t\t\t"
    #     output_line += "\t"+ppi["interAC_propagation"]

    output_line += "\tbinary-interaction"
    #output_line += "\tannotation_id_whatodo"+str(ppi_number)
    output_line += "\tENYO"
    output_line += "\tcurated"
    output_line += "\tPROTEIN"
    output_line += "\tNX_"+ppi["interactor1_accession"]
    output_line += "\tECO:0000353"
    output_line += "\tGOLD"
    #output_line += "\traw_statement_id_whatodo"+str(ppi_number)
    output_line += "\tPubMed"
    output_line += "\tpublication"
    output_line += "\tENYO"
    #output_line += "\tstatement_id_whatodo"
    #output_line += "\ttarget_isoform_howtoknow"
    #output_line += "\t<null>" * 30
    ppi_list.append(output_line)


##Function to store the interaction mapping information
def store_interaction_mapping(ppi, intmap_list):
    output_line = ppi["pmid"]+"\t"+ppi["psimi_id"]
    output_line += "\tNX_"+ppi["interactor1_accession"]+"\t"+ppi["interactor1_name"]
    output_line += "\tNX_"+ppi["interactor2_accession"]+"\t"+ppi["interactor2_name"]

    if (ppi["interactor1_isoform_accession"] != ""):
        output_line += "\t"+ppi["interactor1_mapping_sequence"]+"\t"+ppi["interactor1_isoform_accession"]
        output_line += "\t"+ppi["interactor1_occurrence_start"]+"\t"+ppi["interactor1_occurrence_stop"]+"\t"+ppi["interactor1_occurrence_identity"]
        
    elif (ppi["interactor2_isoform_accession"] != ""):
        output_line += "\t"+ppi["interactor2_mapping_sequence"]+"\t"+ppi["interactor2_isoform_accession"]
        output_line += "\t"+ppi["interactor2_occurrence_start"]+"\t"+ppi["interactor2_occurrence_stop"]+"\t"+ppi["interactor2_occurrence_identity"]

    # if (ppi["ppi_type"] == "vh"):
    #     output_line += "\t"+ppi["Uniprot_viral_species"]+"\t"+ppi["FTId"]+"\t"+ppi["Chain_name"]
    #     output_line += "\t"+ppi["interAC_propagation"]
    # elif (ppi["ppi_type"] == "hh"):
    #     output_line += "\t\t\t"
    #     output_line += "\t"+ppi["interAC_propagation"]

    intmap_list.append(output_line)


##Function to store every interaction mapping information
def integration_formatting(propagated_ppi_file, ppi_output, interaction_mapping_output):
    logging.info("Starts integration formatting")
    ppi_result = open(ppi_output,"w")
    intmap_result = open(interaction_mapping_output, "w")

    ##write the header in PPi output file
    ppi_header = "reference_accession\tpsimi_id"
    ppi_header += "\tnextprot_accession\tgene_name" #interactant1
    ppi_header += "\tbiological_object_accession\tbiological_object_name" #interactant2
    #ppi_header += "\tUniprot_viral_species\tFTId\tChain_name"
    #ppi_header += "\tinterAC_propagation"

    #extra values
    ppi_header += "\tannotation_category"
    #ppi_header += "\tannotation_id"
    ppi_header += "\tassigment_method"
    ppi_header += "\tassigned_by"
    ppi_header += "\tbiological_object_type"
    ppi_header += "\tentry_accession"
    ppi_header += "\tevidence_code"
    ppi_header += "\tevidence_quality"
    #ppi_header += "\traw_statement_id"
    ppi_header += "\treference_database"
    ppi_header += "\tresource_type"
    ppi_header += "\tsource"
    #ppi_header += "\tstatement_id"
    #ppi_header += "\ttarget_isoforms"
    ppi_header += "\tannot_source_accession"

    # #null
    # ppi_header += "\tannotation_name\tannotation_object_species\tannotation_subject_species"
    # ppi_header += "\tannot_cv_term_accession\tannot_cv_term_name\tannot_cv_term_terminology\tannot_description"
    # ppi_header += "\tbiological_object_database\tdebug_info\tevidence_intensity\tevidence_note\tevidence_properties\tevidence_statement_ref"
    # ppi_header += "\texp_context_eco_detect_method\texp_context_eco_iss\exp_context_eco_mutation"
    # ppi_header += "\textra_fields\tisoform_canonical\tis_negative\tlocation_begin\tlocation_begin_master\tlocation_end\tlocation_end_master"
    # ppi_header += "\tobject_annotation_ids\tobject_annot_entry_unames\tobject_annot_iso_unames\tobject_statement_ids"
    # ppi_header += "\tsubject_annotation_ids\tsubject_statement_ids\tvariant_original_amino_acid\tvariant_variation_amino_acid"

    ppi_result.write(ppi_header.upper()+"\n")

    ##write the header in interaction mapping output file
    intmap_header = "pmid\tpsimi_id"
    intmap_header += "\tinteractor1_accession\tinteractor1_name"
    intmap_header += "\tinteractor2_accession\tinteractor2_name"
    intmap_header += "\tinteractor_mapping_sequence\tinteractor_isoform_accession"
    intmap_header += "\tinteractor_occurrence_start\tinteractor_occurrence_stop\tinteractor_occurrence_identity" 
    # intmap_header += "\tUniprot_viral_species\tFTId\tChain_name"
    # intmap_header += "\tinterAC_propagation"
    intmap_result.write(intmap_header.upper()+"\n")
    
    ppi_list=[]
    intmap_list=[]

    with open(propagated_ppi_file) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        ##fill lists of ppi and interaction mapping
        for ppi in reader:
            store_ppi(ppi, ppi_list)
            if ( (ppi["interactor1_isoform_accession"] != "") or (ppi["interactor2_isoform_accession"] != "") ):
                store_interaction_mapping(ppi, intmap_list)

        ##write
        ppi_list = list(set(ppi_list)) ##save one copy of the duplicated PPi
        cpt = 0
        for uniq_ppi in ppi_list:
            cpt += 1
            ppi_result.write(uniq_ppi+"\t"+str(cpt)+"\n")

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
