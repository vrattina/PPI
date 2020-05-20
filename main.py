#AUTHOR: RATTINA Vimel - 2019/11/28 - SIB & Enyo Pharma

import sys #I/O files
import configparser
import os
import shutil
from subprocess import call #integration_formatting is in python3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("enyo.config", help="Configuration file with ENYO settings to fill")
parser.add_argument("--step", nargs='+', help="[all, flatten, download, canonical, update, classify, select, propagation, symmetry, integration, json]", choices=["all", "flatten", "download", "canonical", "update", "classify", "select", "propagation", "symmetry", "integration", "json"])
args = parser.parse_args()


#if len(sys.argv) < 2:
#    sys.exit('Usage: %s <configuration_file> <step>\n<configuration_file>:\n<step>:' % sys.argv[0])

# from flatten_enyo_json import flatten_enyo_json
# from download_ac_in_uniprot_versions import download_ac_in_uniprot_versions
# from find_canonical_isoform import find_canonical_isoform
# from update_ac_uniprot_information import update_ac_uniprot_information #dependency varsplic.pl & swissknife
# from classify_ppi import classify_ppi
# from select_one_isoform import select_one_isoform
# from inter_ac_propagation import inter_ac_propagation #dependency human_identical_proteins.sparql & neXtProt API
# from symmetric_interaction import symmetric_interaction
# ##from integration_formatting import integration_formatting #cannot specify option
# from tsv_to_json import tsv_to_json

##The aim of this script is to:
##-Transform a txt file with 1 to n json VH and or HH-PPi into human readable flat file,
##-Download the viral and human proteins tab (history) and the AC.txt of the 2 UniProt version of interest files,
##-Compare the gene name / sequence between the 2 AC.txt UniProt version to notify if there is any changes,
##-Add information for the viral proteins (viral species and chain name + FTId if corresponding viral positions).

## OUTPUTS ##
##-Two files to be integrated in neXtProt (one with PPi and one with interaction mapping),
##-A file to be re-curated by neXtProt and ENYO pharma before to be integrated.

def main(configuration_file):

    Config = configparser.ConfigParser()
    Config.read(configuration_file)
    
    ppi_in_json = Config.get('ENYO', 'ppi_in_json')
    uniprot_version_enyo = Config.get('ENYO', 'uniprot_version_enyo')
    uniprot_version_nextprot = Config.get('ENYO', 'uniprot_version_nextprot')
    proxy_folder = Config.get('ENYO', 'proxy_folder')
    tmp_folder = Config.get('ENYO', 'tmp_folder')
    output_folder = Config.get('ENYO', 'output_folder')
    
    if not os.path.exists(proxy_folder):
        os.makedirs(proxy_folder)

    if not os.path.exists(tmp_folder):
        os.makedirs(tmp_folder)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    save_input = shutil.copy(ppi_in_json, proxy_folder)

    flatten_ppi_tsv = tmp_folder + "/enyo_flatten.tsv"
    if "flatten" in args.step or args.step == "all":
        exit_code = call("python3.6 flatten_enyo_json.py " + ppi_in_json + " " + flatten_ppi_tsv, shell=True)
        
    if "download" in args.step or args.step == "all":
        exit_code = call("python3.6 download_ac_in_uniprot_versions.py " + flatten_ppi_tsv + " " + uniprot_version_enyo + " " + proxy_folder + " True", shell=True)
        exit_code = call("python3.6 download_ac_in_uniprot_versions.py " + flatten_ppi_tsv + " " + uniprot_version_nextprot + " " + proxy_folder + " False", shell=True)

    enyo_folder = proxy_folder + "/" + uniprot_version_enyo
    flatten_ppi_specified_iso_tsv = tmp_folder + "/enyo_flatten_specified_iso.tsv"
    if "canonical" in args.step or args.step == "all":
        exit_code = call("python3.6 find_canonical_isoform.py " + flatten_ppi_tsv + " " + enyo_folder + " " + flatten_ppi_specified_iso_tsv, shell=True)

    VP_dict = tmp_folder + "/VPdict.dict"
    HP_dict = tmp_folder + "/HPdict.dict"
    if "update" in args.step or args.step == "all":
        exit_code = call("python3.6 update_ac_uniprot_information.py " + flatten_ppi_specified_iso_tsv + " " + proxy_folder + " " + uniprot_version_enyo + " " + uniprot_version_nextprot + " " + VP_dict + " " + HP_dict, shell=True)

    raw_integrate_file = tmp_folder + "/integrate_raw.tsv"
    outputfile_curate = output_folder + "/ppi_to_cure.tsv"
    if "classify" in args.step or args.step == "all":
        exit_code = call("python3.6 classify_ppi.py " + flatten_ppi_specified_iso_tsv + " " + VP_dict + " " + HP_dict + " " + raw_integrate_file + " " + outputfile_curate, shell=True)

    one_isoform_file = tmp_folder + "/integrate_one_iso.tsv"
    if "select" in args.step or args.step == "all":
        exit_code = call("python3.6 select_one_isoform.py " + raw_integrate_file + " " + one_isoform_file, shell=True)

    propagated_file = tmp_folder + "/integrate_one_iso_propagated.tsv"
    if "propagation" in args.step or args.step == "all":
        exit_code = call("python3.6 inter_ac_propagation.py " + one_isoform_file + " " + propagated_file, shell=True)

    symmetric_file = tmp_folder + "/integrate_one_iso_propagated_symmetric.tsv"
    if "symmetry" in args.step or args.step == "all":
        exit_code = call("python3.6 symmetric_interaction.py " + propagated_file + " " + symmetric_file, shell=True)

    ppi_output = tmp_folder + "/integrate_binaryinteraction.tsv"
    intmap_output = tmp_folder + "/integrate_intmap.tsv"
    if "integration" in args.step or args.step == "all":
        exit_code = call("python3.6 integration_formatting.py " + symmetric_file + " " + ppi_output + " " + intmap_output + " psimi_unmerged", shell=True)

    ppi_json = output_folder + "/integrate_binaryinteraction.json"
    intmap_json = output_folder + "/integrate_intmap.json"
    if "json" in args.step or args.step == "all":
        exit_code = call("python3.6 tsv_to_json.py " + ppi_output + " " + ppi_json, shell=True)
        exit_code = call("python3.6 tsv_to_json.py " + intmap_output + " " + intmap_json, shell=True)

if __name__ == "__main__":
    main(sys.argv[1])

#python3.6 main.py /home/vrattina/Documents/HH_PPi/enyo.config --step all
