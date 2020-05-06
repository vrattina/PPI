#!/usr/bin/python
#AUTHOR: RATTINA Vimel - 2020/04/29 - SIB & Enyo Pharma

import sys #I/O files
import json #json loader
import logging #log file
import csv
import re

logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', filename='main.log', level=logging.DEBUG)

##The aim of this script is to focus on interaction mapping, found the isoform number of the canonical isoform (not repertoried by ENYO) to update it

## INPUTS: ##
##-The flatten PPi file from ENYO (output of flatten_enyo_json.py),
##-ENYO YYYY_MM UniProt version folder pathway where the AC.txt are stored

## OUTPUTS: ##
##-The flatten PPi file from ENYO with the canonical isoform specification


###################
# FUNCTIONS #######
###################

##function seeking for the canonical isoform in ENYO YYYY_MM UniProt version
def make_dictionnary_canonical_isoform(flatten_ppi_tsv, YYYY_MM_enyo_folder):

    #First create a list of unique AC of isoform in interaction mapping
    isoform_list = []

    with open(flatten_ppi_tsv) as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")        
        for ppi in reader:
            int1iso = ppi["interactor1_isoform_accession"]
            int2iso = ppi["interactor2_isoform_accession"]
    
            #if no isoform information, the canonical form was taken, seeking for it
            if (int1iso != "") and ("-" not in int1iso):
                isoform_list.append(int1iso)

            if (int2iso != "") and ("-" not in int2iso):
                isoform_list.append(int2iso)
    isoform_list = list(set(isoform_list))
    
    #Then look in ENYO version AC.txt the canonical sequence displayed
    canonical_dict={}

    for iso in isoform_list:
        iso_path =YYYY_MM_enyo_folder+"/"+iso+".txt" 
        with open(iso_path) as iso_enyo:
            full_txt =  iso_enyo.read()
            canonical_pattern = "IsoId=(" + re.escape(iso) + "-\d+); Sequence=Displayed;"
            canonical_AC = re.search(canonical_pattern, full_txt, re.DOTALL)
            
            #if the canonical information is diplayed took it
            if canonical_AC is not None:
                canonical_dict[iso] = canonical_AC.group(1)

            #if not displayed, there is no alternative isoform or it is -1 so not specified
            if canonical_AC is None:
                canonical_dict[iso] = iso+"-1"

    return canonical_dict


def find_canonical_isoform(flatten_ppi_tsv, YYYY_MM_enyo_folder, outputfile):
    
    logging.info("Seeking for the canonical isoform AC (interaction mapping) in ENYO UniProt version")

    canonical_dict = make_dictionnary_canonical_isoform(flatten_ppi_tsv, YYYY_MM_enyo_folder)
    #print str(canonical_dict)

    logging.info("Write the canonical isoform AC in interaction mapping when needed")

    outfile = open(outputfile, "w")

    with open(flatten_ppi_tsv) as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")

        
        header_list = reader.fieldnames
        header = "\t".join(str(elt) for elt in header_list)
        header += "\n"

        outfile.write(header)

        for ppi in reader:
         
            ppitype = ppi["ppi_type"]
            stableid = ppi["stable_id"]
            pmid = ppi["pmid"]
            psimi = ppi["psimi_id"]
                    
            int1ac = ppi["interactor1_accession"]
            int1name = ppi["interactor1_name"]
            int1start = ppi["interactor1_start"]
            int1stop = ppi["interactor1_stop"]

            int2ac = ppi["interactor2_accession"]
            int2name = ppi["interactor2_name"]
            int2start = ppi["interactor2_start"]
            int2stop = ppi["interactor2_stop"]
    
            int1mapseq = ppi["interactor1_mapping_sequence"]
            int1iso = ppi["interactor1_isoform_accession"]
            int1occstart = ppi["interactor1_occurrence_start"]
            int1occstop = ppi["interactor1_occurrence_stop"]
            int1occidentity = ppi["interactor1_occurrence_identity"]
    
            int2mapseq = ppi["interactor2_mapping_sequence"]
            int2iso = ppi["interactor2_isoform_accession"]
            int2occstart = ppi["interactor2_occurrence_start"]
            int2occstop = ppi["interactor2_occurrence_stop"]
            int2occidentity = ppi["interactor2_occurrence_identity"]

            #if isoform found in the canonical dictionnary replace it by the value
            if int1iso in canonical_dict:
                int1iso = canonical_dict[int1iso]

            if int2iso in canonical_dict:
                int2iso = canonical_dict[int2iso]

            line = ppitype+"\t"+stableid+"\t"+pmid+"\t"+psimi
            line += "\t"+int1ac+"\t"+int1name+"\t"+int1start+"\t"+int1stop
            line += "\t"+int2ac+"\t"+int2name+"\t"+int2start+"\t"+int2stop
            line += "\t"+int1mapseq+"\t"+int1iso
            line += "\t"+int1occstart+"\t"+int1occstop+"\t"+int1occidentity
            line += "\t"+int2mapseq+"\t"+int2iso
            line += "\t"+int2occstart+"\t"+int2occstop+"\t"+int2occidentity+"\n"

            outfile.write(line)

    outfile.close()

    logging.info("End of canonical isoform AC writting")

####################
# MAIN #############
####################

if len(sys.argv) < 4:
    sys.exit('Usage: %s <flat_file> <MM_YYYY_ENYO_folder> <output_file>\n<flat_file>: A flatten tsv file where a new line is created when one information in the json changes i.e. when the interaction mapping is available because it could be found many times in the same isoform [position changes] and/or in many isoforms [protein AC changes] with more or less the same identity\n<MM_YYYY_ENYO_folder>: folder where AC.txt of ENYO are stored\n<output_file>:A flatten tsv file with canonical isoform specified' % sys.argv[0])

if __name__ == "__main__":
    find_canonical_isoform(sys.argv[1], sys.argv[2], sys.argv[3])
#./find_canonical_isoform.py vinland-hh-2020-01-22.tsv /home/vrattina/Documents/HH_PPi/hh-2020-04-22_results/2019_01/ iso_vinland-hh-2020-01-22.tsv
