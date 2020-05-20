#AUTHOR: RATTINA Vimel - 2020/02/07 - SIB & Enyo Pharma

import sys #I/O files
import os, errno #create folder
import csv #parse tsv
import logging #log file

logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', filename='main.log', level=logging.DEBUG)

###The aim of this script is to permit to have a symmetric PPi file

## INPUTS: ##
##-The third-version file (3/4) with propagated PPi

## OUTPUTS: ##
##-An output file with propagated and symmetric PPi, the fourth-version of PPi intermediar file

###################
# FUNCTIONS #######
###################

##Function writting the symmetric PPi
def symmetric_interaction(one_isoform_propagated_file, propagated_symmetric_ppi_output):
    logging.info("Starts writting symmetric PPi")
    with open(one_isoform_propagated_file) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')

        output = open(propagated_symmetric_ppi_output, "w")

        header_list = reader.fieldnames
        header = "\t".join(str(elt) for elt in header_list)
        header += "\n"
        output.write(header)

        symmetric_ppi_list=[]
        for ppi in reader:

            description = ppi["ppi_type"]+"\t"+ppi["pmid"]+"\t"+ppi["psimi_id"]

            interaction1_part1 = "\t"+ppi["interactor1_accession"]+"\t"+ppi["interactor1_name"]
            interaction2_part1 = "\t"+ppi["interactor2_accession"]+"\t"+ppi["interactor2_name"]

            interaction1_part2 = "\t"+ppi["interactor1_mapping_sequence"]
            interaction1_part2 += "\t"+ppi["interactor1_isoform_accession"]
            interaction1_part2 += "\t"+ppi["interactor1_occurrence_start"]
            interaction1_part2 += "\t"+ppi["interactor1_occurrence_stop"]
            interaction1_part2 += "\t"+ppi["interactor1_occurrence_identity"]
            
            interaction2_part2 = "\t"+ppi["interactor2_mapping_sequence"]
            interaction2_part2 += "\t"+ppi["interactor2_isoform_accession"]
            interaction2_part2 += "\t"+ppi["interactor2_occurrence_start"]
            interaction2_part2 += "\t"+ppi["interactor2_occurrence_stop"]
            interaction2_part2 += "\t"+ppi["interactor2_occurrence_identity"]

            viral_info = "\t"+ppi["Uniprot_viral_species"]+"\t"+ppi["FTId"]+"\t"+ppi["Chain_name"]
            propagation_info = "\t"+ppi["interAC_propagation"]+"\n"

            #interactant1 before interactant2
            ppi_symmetry1 = description
            ppi_symmetry1 += interaction1_part1+interaction2_part1
            ppi_symmetry1 += interaction1_part2+interaction2_part2
            ppi_symmetry1 += viral_info+propagation_info

            #interactant2 before interactant1
            ppi_symmetry2 = description
            ppi_symmetry2 += interaction2_part1+interaction1_part1
            ppi_symmetry2 += interaction2_part2+interaction1_part2
            ppi_symmetry2 += viral_info+propagation_info

            symmetric_ppi_list.append(ppi_symmetry1)

            if (ppi_symmetry1 != ppi_symmetry2): 
                symmetric_ppi_list.append(ppi_symmetry2)

        #for some homodimer enyo already done symmetric ppi
        symmetric_ppi_list = list(set(symmetric_ppi_list))
        for i in symmetric_ppi_list:
            output.write(str(i))

        logging.info("Writting symmetric PPi ends")

####################
# MAIN #############
####################

if len(sys.argv) < 3:
    sys.exit('Usage: %s <propagated_ppi_file> <propagated_symmetric_output_file> \n<propagated_ppi_file>: third version of the file to be integrated in neXtProt\n<propagated_symmetric_output_file>: fourth version of the file to be integrated in neXtProt (PPi are writting in symmetric)' % sys.argv[0])

if __name__ == "__main__":
    symmetric_interaction(sys.argv[1], sys.argv[2])

#./symmetric_interaction.py one_iso_integre_vinland-hh-2020-01-22.tsv ppi_integre_vinland-hh-2020-01-22.tsv intmap_integre_vinland-hh-2020-01-22.tsv
