#AUTHOR: RATTINA Vimel - 2020/02/12 - SIB & Enyo Pharma

import sys #I/O files
import os, errno #create folder
import csv #parse tsv
import logging #log file
import json
import csv

logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', filename='main.log', level=logging.DEBUG)

###Convert a tsv file to json format

## INPUTS: ##
##-A tsv file
## OUTPUTS: ##
##-A json file

###################
# FUNCTIONS #######
###################

##Function

##Function writting the output file with the propagated PPi
def tsv_to_json(tsv_file, json_output):
    logging.info("Starts "+tsv_file+" conversion into json")
    
    json_file = open(json_output, "w")

    with open(tsv_file) as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        data = list(reader)
    
        json_data = json.dumps(data)
        json_file.write(json_data)

    json_file.close()
    logging.info("Conversion into json ends")

####################
# MAIN #############
####################

if len(sys.argv) < 3:
    sys.exit('Usage: %s <tsv_file> <json_output>\n<tsv_file>: a tsv file\n<json_output>: a json file' % sys.argv[0])

if __name__ == "__main__":
    tsv_to_json(sys.argv[1], sys.argv[2])

#./tsv_to_json.py ppi_ntintegre ppi_ntintegre.json
