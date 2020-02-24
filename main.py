#!/usr/bin/python
#AUTHOR: RATTINA Vimel - 2019/11/28 - SIB & Enyo Pharma

import sys #I/O files

if len(sys.argv) < 8:
    sys.exit('Usage: %s <ppi_in_json> <flatten_ppi_tsv> <YYYY_MM_enyo> <YYYY_MM_nextprot> <outputfolder> <outputfile_integrate> <outputfile_curate>\n<ppi_in_json>:\n<flatten_ppi_tsv>:\n<YYYY_MM_enyo>:\n<YYYY_MM_nextprot>:\n<outputfolder>\n<outputfile_integrate>:\n<outputfile_curate>:' % sys.argv[0])

from download_ac_in_uniprot_versions import download_ac_in_uniprot_versions
from flatten_enyo_json import flatten_enyo_json
from update_ac_uniprot_information import update_ac_uniprot_information #dependency varsplic.pl & swissknife
from classify_ppi import classify_ppi
from select_one_isoform import select_one_isoform
from inter_ac_propagation import inter_ac_propagation #dependency human_identical_proteins.sparql & neXtProt API
from integration_formatting import integration_formatting
from tsv_to_json import tsv_to_json

##The aim of this script is to:
##-Transform a txt file with 1 to n json VH and or HH-PPi into human readable flat file,
##-Download the viral and human proteins tab (history) and the AC.txt of the 2 UniProt version of interest files,
##-Compare the gene name / sequence between the 2 AC.txt UniProt version to notify if there is any changes,
##-Add information for the viral proteins (viral species and chain name + FTId if corresponding viral positions).

## OUTPUTS ##
##-Two files to be integrated in neXtProt (one with PPi and one with interaction mapping),
##-A file to be re-curated by neXtProt and ENYO pharma before to be integrated.

def main(ppi_in_json, flatten_ppi_tsv, YYYY_MM_enyo, YYYY_MM_nextprot, outputfolder, outputfile_integrate, outputfile_curate):

    #flatten_enyo_json(ppi_in_json, flatten_ppi_tsv)
    
    #download_ac_in_uniprot_versions(flatten_ppi_tsv, YYYY_MM_enyo, outputfolder, "True")
    #download_ac_in_uniprot_versions(flatten_ppi_tsv, YYYY_MM_nextprot, outputfolder, "False") 

    VP_dict = "VPdict_"+flatten_ppi_tsv.replace(".tsv",".dict")
    HP_dict = "HPdict_"+flatten_ppi_tsv.replace(".tsv",".dict")
    #update_ac_uniprot_information(flatten_ppi_tsv, outputfolder, YYYY_MM_enyo, YYYY_MM_nextprot, VP_dict, HP_dict)
    
    raw_integrate_file = outputfile_integrate+"_raw"
    classify_ppi(flatten_ppi_tsv, VP_dict, HP_dict, raw_integrate_file, outputfile_curate)

    one_isoform_file = outputfile_integrate+"_one_iso"
    select_one_isoform(raw_integrate_file, one_isoform_file)

    propagated_file = one_isoform_file+"_propagated"
    inter_ac_propagation(one_isoform_file, propagated_file)

    ppi_output = "ppi_"+outputfile_integrate
    intmap_output = "intmap_"+outputfile_integrate
    integration_formatting(propagated_file, ppi_output, intmap_output)

    ppi_json = ppi_output+".json"
    intmap_json = intmap_output+".json"
    tsv_to_json(ppi_output, ppi_json)
    tsv_to_json(intmap_output, intmap_json)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])

#./main.py vinland-vh-hh-ppi-sib vinland-vh-hh-ppi-sib.tsv 2019_01 2019_11 vh_hh_results vh_hh_integre.tsv vh_hh_cure.tsv &> vh_hh.log
#./main.py vinland-hh-2020-01-22 vinland-hh-2020-01-22.tsv 2019_01 2019_11 hh-2020-01-22_results/ integre_vinland-hh-2020-01-22.tsv cure_vinland-hh-2020-01-22.tsv $> hh.log 

#To improve:
#if AC replaced by a new AC due to trEmbl to SwissProt --> notify it

