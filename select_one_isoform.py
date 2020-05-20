#AUTHOR: RATTINA Vimel - 2020/02/07 - SIB & Enyo Pharma

import sys #I/O files
import os, errno #create folder
import csv #parse tsv
import logging #log file

logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', filename='main.log', level=logging.DEBUG)

###The aim of this script is to only select the interesting information from the file to be integrated, and it select only one isoform by interactant couple+psi-mi+pmid and interaction mapping.

## INPUTS: ##
##-The first-version file (1/4) to be integrated by neXtProt with no gene name or sequence changes between 2 versions of uniprot i.e. the flatten PPi with in addition viral species, chain name and FTId

## OUTPUTS: ##
##-An output file (2/4) with the interesting information to be integrated in neXtProt and with only one isoform (isoform propagation will be handle by neXtProt inner script)

###################
# FUNCTIONS #######
###################

##Function writting in the output with the requirred information
def output_writter(ppi, output_file):
    output_line = ppi["ppi_type"]+"\t"+ppi["pmid"]+"\t"+ppi["psimi_id"]
    output_line += "\t"+ppi["interactor1_accession"]+"\t"+ppi["interactor1_name"]
    output_line += "\t"+ppi["interactor2_accession"]+"\t"+ppi["interactor2_name"]
    output_line += "\t"+ppi["interactor1_mapping_sequence"]+"\t"+ppi["interactor1_isoform_accession"]
    output_line += "\t"+ppi["interactor1_occurrence_start"]+"\t"+ppi["interactor1_occurrence_stop"]+"\t"+ppi["interactor1_occurrence_identity"]
    output_line += "\t"+ppi["interactor2_mapping_sequence"]+"\t"+ppi["interactor2_isoform_accession"]
    output_line += "\t"+ppi["interactor2_occurrence_start"]+"\t"+ppi["interactor2_occurrence_stop"]+"\t"+ppi["interactor2_occurrence_identity"]
    output_line += "\t"+ppi["Uniprot_viral_species"]+"\t"+ppi["FTId"]+"\t"+ppi["Chain_name"]+"\n"
    output_file.write(output_line)

##Function writting the output file with only one occurrence to be propagated on other isoforms in neXtProt
def select_one_isoform(raw_integrate_file, one_isoform_output):
    logging.info("Starts selecting one occurrence")
    output_file = open(one_isoform_output,"w")

    ##write the header in both output files
    header = "ppi_type\tpmid\tpsimi_id"

    header += "\tinteractor1_accession\tinteractor1_name"
    header += "\tinteractor2_accession\tinteractor2_name"

    header += "\tinteractor1_mapping_sequence\tinteractor1_isoform_accession"
    header += "\tinteractor1_occurrence_start\tinteractor1_occurrence_stop\tinteractor1_occurrence_identity" 

    header += "\tinteractor2_mapping_sequence\tinteractor2_isoform_accession"
    header += "\tinteractor2_occurrence_start\tinteractor2_occurrence_stop\tinteractor2_occurrence_identity" 

    header +="\tUniprot_viral_species\tFTId\tChain_name"

    output_file.write(header+"\n")

    with open(raw_integrate_file) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')

        ##write the first line
        first_ppi = next(reader)

        output_writter(first_ppi, output_file)
        
        ##variable initialization
        prev_intmap_interactor1 = first_ppi["interactor1_mapping_sequence"]
        prev_iso_interactor1 = first_ppi["interactor1_isoform_accession"]
        prev_occstart_interactor1 = first_ppi["interactor1_occurrence_start"]
        prev_occstop_interactor1 = first_ppi["interactor1_occurrence_stop"]

        prev_intmap_interactor2 = first_ppi["interactor2_mapping_sequence"]
        prev_iso_interactor2 = first_ppi["interactor2_isoform_accession"]
        prev_occstart_interactor2 = first_ppi["interactor2_occurrence_start"]
        prev_occstop_interactor2 = first_ppi["interactor2_occurrence_stop"]
    
        prev_stableid = first_ppi["stable_id"]
        prev_interactor1 = first_ppi["interactor1_accession"]
        prev_interactor2 = first_ppi["interactor2_accession"] 
        print(prev_intmap_interactor1+"\t"+prev_iso_interactor1+"\t"+prev_occstart_interactor1+"\t"+prev_occstop_interactor1)
        for ppi in reader:

            ##if there is no interaction mapping or if the stable id is different (so at least one interactor, pmid or psi-mi)
            ##can update the previous values and write the result
            if ( (ppi["interactor1_isoform_accession"] == "") and (ppi["interactor2_isoform_accession"] == "") ) or ( ppi["stable_id"] != prev_stableid ):

                output_writter(ppi, output_file)

                prev_stableid = ppi["stable_id"]
                prev_interactor1 = ppi["interactor1_accession"]
                prev_interactor2 = ppi["interactor2_accession"]

                prev_intmap_interactor1 = ppi["interactor1_mapping_sequence"]
                prev_iso_interactor1 = ppi["interactor1_isoform_accession"]
                prev_occstart_interactor1 = ppi["interactor1_occurrence_start"]
                prev_occstop_interactor1 = ppi["interactor1_occurrence_stop"]
        
                prev_intmap_interactor2 = ppi["interactor2_mapping_sequence"]
                prev_iso_interactor2 = ppi["interactor2_isoform_accession"]
                prev_occstart_interactor2 = ppi["interactor2_occurrence_start"]
                prev_occstop_interactor2 = ppi["interactor2_occurrence_stop"]

                continue ##continue permit to ignore the remaining code and go to the next iterative loop

            ##if there is an interaction mapping in one interactant, the mapping sequence is identical but not found in the same position (different annotation)
            ##can update the previous values and write the result
            if ( (prev_iso_interactor1 == ppi["interactor1_isoform_accession"]) \
                 and (ppi["interactor1_isoform_accession"] != "") \
                 and ( ppi["interactor1_mapping_sequence"] == prev_intmap_interactor1 ) \
                 and ( ppi["interactor1_occurrence_start"] != prev_occstart_interactor1 ) \
                 and ( ppi["interactor1_occurrence_stop"] != prev_occstop_interactor1 ) ) \
                or \
                ( (prev_iso_interactor2 == ppi["interactor2_isoform_accession"]) \
                  and (ppi["interactor2_isoform_accession"] != "") \
                  and ( ppi["interactor2_mapping_sequence"] == prev_intmap_interactor2 ) \
                  and ( ppi["interactor2_occurrence_start"] != prev_occstart_interactor2 ) \
                  and ( ppi["interactor2_occurrence_stop"] != prev_occstop_interactor2 ) ) :

                output_writter(ppi, output_file)

                prev_stableid = ppi["stable_id"]
                prev_interactor1 = ppi["interactor1_accession"]
                prev_interactor2 = ppi["interactor2_accession"]

                prev_intmap_interactor1 = ppi["interactor1_mapping_sequence"]
                prev_iso_interactor1 = ppi["interactor1_isoform_accession"]
                prev_occstart_interactor1 = ppi["interactor1_occurrence_start"]
                prev_occstop_interactor1 = ppi["interactor1_occurrence_stop"]
        
                prev_intmap_interactor2 = ppi["interactor2_mapping_sequence"]
                prev_iso_interactor2 = ppi["interactor2_isoform_accession"]
                prev_occstart_interactor2 = ppi["interactor2_occurrence_start"]
                prev_occstop_interactor2 = ppi["interactor2_occurrence_stop"]

                continue ##continue permit to ignore the remaining code and go to the next iterative loop
                
            ##knowing that the interaction mapping are ordered by canonical to alternative isoform, if the previous is the canonical do not write the alternative isoform
            ##not in permit to select only the canonical isoform or the first alternative isoform
            ##if we are here in the code, the both interaction mapping are not empty, check if interactor isoform 1 is not empty otherwise write nothing
            prev_AC1 = ""
            if prev_iso_interactor1 != "":
                prev_AC1 = prev_iso_interactor1.split("-")[0]

            if ( prev_AC1 not in ppi["interactor1_isoform_accession"] ) \
               and ( ppi["interactor1_isoform_accession"] != "" ) \
               and ( ppi["interactor1_mapping_sequence"] == prev_intmap_interactor1 ):

                prev_stableid = ppi["stable_id"]
                prev_interactor1 = ppi["interactor1_accession"]
                prev_interactor2 = ppi["interactor2_accession"]

            ##if isoform/mapping not empty and interaction mapping sequence is different from the previous one, write
            elif ( ppi["interactor1_mapping_sequence"] != prev_intmap_interactor1 ):

                output_writter(ppi, output_file)

                prev_intmap_interactor1 = ppi["interactor1_mapping_sequence"]
                prev_iso_interactor1 = ppi["interactor1_isoform_accession"]
                prev_occstart_interactor1 = ppi["interactor1_occurrence_start"]
                prev_occstop_interactor1 = ppi["interactor1_occurrence_stop"]
            
                prev_intmap_interactor2 = ppi["interactor2_mapping_sequence"]
                prev_iso_interactor2 = ppi["interactor2_isoform_accession"]
                prev_occstart_interactor2 = ppi["interactor2_occurrence_start"]
                prev_occstop_interactor2 = ppi["interactor2_occurrence_stop"]

                continue

            ##same for interactor2
            prev_AC2 = ""
            if prev_iso_interactor1 != "":
                prev_AC2 = prev_iso_interactor1.split("-")[0]

            if ( prev_AC2 not in ppi["interactor2_isoform_accession"] ) \
               and ( ppi["interactor2_isoform_accession"] != "" ) \
               and ( ppi["interactor2_mapping_sequence"] == prev_intmap_interactor2 ):
                prev_stableid = ppi["stable_id"]
                prev_interactor1 = ppi["interactor1_accession"]
                prev_interactor2 = ppi["interactor2_accession"]  

            elif( ppi["interactor2_mapping_sequence"] != prev_intmap_interactor2 ):

                output_writter(ppi, output_file)
                prev_intmap_interactor1 = ppi["interactor1_mapping_sequence"]
                prev_iso_interactor1 = ppi["interactor1_isoform_accession"]
                prev_occstart_interactor1 = ppi["interactor1_occurrence_start"]
                prev_occstop_interactor1 = ppi["interactor1_occurrence_stop"]

            
                prev_intmap_interactor2 = ppi["interactor2_mapping_sequence"]
                prev_iso_interactor2 = ppi["interactor2_isoform_accession"]
                prev_occstart_interactor2 = ppi["interactor2_occurrence_start"]
                prev_occstop_interactor2 = ppi["interactor2_occurrence_stop"]

                continue
    
    output_file.close()
    logging.info("Selecting one occurrence ends")

####################
# MAIN #############
####################

if len(sys.argv) < 3:
    sys.exit('Usage: %s <raw_integrate_file> <one_isoform_output>\n<raw_integrate_file>: first version of the file to be integrated in neXtProt\n<one_isoform_output>: second version of the file to be integrated in neXtProt (only one isoform is selected for each PPi and interaction mapping)\n' % sys.argv[0])

if __name__ == "__main__":
    select_one_isoform(sys.argv[1], sys.argv[2])

#./select_one_isoform.py raw_integre_vinland-hh-2020-01-22.tsv one_iso_integre_vinland-hh-2020-01-22.tsv
