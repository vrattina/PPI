#AUTHOR: RATTINA Vimel - 2020/01/31 - SIB & Enyo Pharma

import sys #I/O files
import json #json loader
import logging #log file

logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', filename='main.log', level=logging.DEBUG)

##The aim of this script is to transform a PPi file with 1 to n json PPi (provided by ENYO pharma) into human readable flat file that could be easily computed

###################
# FUNCTIONS #######
###################

### Function printing the occurrence if only one interactant has interaction mapping information (e.g. every VH-PPi)
def flatten_one_mapping(ppi, interactor, common_core, outfile):
    empty_mapping="\t\t\t\t\t"
    for mapping in ppi[interactor]["mapping"]:
        #print "MAPPING\n"+str(mapping)+"\n"
        
        if (type(mapping) == dict):
            newmapping="\t"+mapping["sequence"]

            for isoform in mapping["isoforms"]:
                newisoform=newmapping+"\t"+str(isoform["accession"])
                #print "ISOFORM\n"+newisoform+"\n" 
                
                for occurrence in isoform["occurrences"]:
                    newoccurrence=newisoform+"\t"+str(occurrence["start"])+"\t"+str(occurrence["stop"])+"\t"+str(occurrence["identity"])#+"\n"

                    ##to have an identical number of columns, add empty mapping for the interactor without mapping
                    if interactor == "interactor1":
                        outfile.write(common_core+newoccurrence+empty_mapping+"\n")

                    elif interactor == "interactor2":
                        outfile.write(common_core+empty_mapping+newoccurrence+"\n")

        else: ##if mapping is not a dict that means enyo json file contains errors
            logging.warning("This description is not treated:\t"+str(ppi))
        

def flatten_enyo_json(jsonfile, outputfile):
    
    logging.info("Converting the ENYO PPi file with several json's into a flat file")


    jfile = open(jsonfile, 'r')
    outfile = open(outputfile, "w")

    ##in ENYO PPi raw data, interactor2 will be the viral protein in vh type PPi
    header = "ppi_type\tstable_id\tpmid\tpsimi_id"
    
    header += "\tinteractor1_accession\tinteractor1_name\tinteractor1_start\tinteractor1_stop"
    header += "\tinteractor2_accession\tinteractor2_name\tinteractor2_start\tinteractor2_stop"
    
    header += "\tinteractor1_mapping_sequence\tinteractor1_isoform_accession"
    header += "\tinteractor1_occurrence_start\tinteractor1_occurrence_stop\tinteractor1_occurrence_identity"
    
    header += "\tinteractor2_mapping_sequence\tinteractor2_isoform_accession"
    header += "\tinteractor2_occurrence_start\tinteractor2_occurrence_stop\tinteractor2_occurrence_identity\n"
    
    outfile.write(header)

    for i in jfile.readlines():

        i = i.rstrip()
        ppi = json.loads(i) ##one json line is a PPi
        #print str(ppi)+"\n"

        ##not every PPi has interaction mapping information, however they all share a minimum of information: a common core
        common_core=ppi["type"]+"\t"+ppi["stable_id"]+"\t"+str(ppi["publication"]["pmid"])+"\t"+str(ppi["method"]["psimi_id"])
        common_core += "\t"+ppi["interactor1"]["protein"]["accession"]+"\t"+ppi["interactor1"]["name"]+"\t"+str(ppi["interactor1"]["start"])+"\t"+str(ppi["interactor1"]["stop"])
        common_core += "\t"+ppi["interactor2"]["protein"]["accession"]+"\t"+ppi["interactor2"]["name"]+"\t"+str(ppi["interactor2"]["start"])+"\t"+str(ppi["interactor2"]["stop"])
        #print "COMMONCORE\n"+common_core+"\n"
        
        empty_mapping="\t\t\t\t\t"
        ##if no mapping, ready to be written in the output file
        if (  (bool(ppi["interactor1"]["mapping"]) == False) and (bool(ppi["interactor2"]["mapping"]) == False)):
            common_core+=empty_mapping+empty_mapping+"\n"
            outfile.write(common_core)

        ##if only one interactor has interaction mapping, use flatten_one_mapping function
        elif (  (bool(ppi["interactor1"]["mapping"]) == True) and (bool(ppi["interactor2"]["mapping"]) == False)):
            flatten_one_mapping(ppi, "interactor1", common_core, outfile)
            
        elif (  (bool(ppi["interactor1"]["mapping"]) == False) and (bool(ppi["interactor2"]["mapping"]) == True)):
            flatten_one_mapping(ppi, "interactor2", common_core, outfile)

        ##if both interactors have interaction mapping, use flatten_two_mapping function            
        elif (  (bool(ppi["interactor1"]["mapping"]) == True) and (bool(ppi["interactor2"]["mapping"]) == True)):
            flatten_one_mapping(ppi, "interactor1", common_core, outfile)
            flatten_one_mapping(ppi, "interactor2", common_core, outfile)
            
    logging.info("ENYO PPi file converted in a flatten file")
    outfile.close()
    jfile.close()

####################
# MAIN #############
####################

if len(sys.argv) < 3:
    sys.exit('Usage: %s <jsons_in_txt_file> <output_file>\n<jsons_in_txt_file>: PPi ENYO output file i.e. a txt file with 1 json for 1 PPi description (P1 interacts with P2 thanks to this psi-mi method in this pmid )\n<output_file>: A flatten tsv file where a new line is created when one information in the json changes i.e. when the interaction mapping is available because it could be found many times in the same isoform [position changes] and/or in many isoforms [protein AC changes] with more or less the same identity' % sys.argv[0])

if __name__ == "__main__":
    flatten_enyo_json(sys.argv[1], sys.argv[2])
#./flatten_enyo_json.py vinland-2019-10-04-sib vinland_flatfile.tsv
