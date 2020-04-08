#!/usr/bin/python3.6
#AUTHOR: RATTINA Vimel - 2020/04/02 - SIB & Enyo Pharma

import sys #I/O files
import os, errno #create folder
import csv #parse tsv
import logging #log file
import argparse

# if len(sys.argv) < 4:
#     sys.exit('Usage: %s <propagated_ppi_file> <ppi_output_file> <interaction_mapping_output_file> <step>\n<propagated_ppi_file>: third version of the file to be integrated in neXtProt\n<ppi_output_file>: PPi file to be integrated in neXtProt\n<interaction_mapping_output_file>: interaction mapping file to be integrated in neXtProt\n<step[psimi_merged, psimi_unmerged]>: psimi_merged to concatenate psimi or psimi_unmerged new psimi new line' % sys.argv[0])

##arguments declaration
parser = argparse.ArgumentParser()
parser.add_argument("propagated_ppi_file", help="third version of the file to be integrated in neXtProt")
parser.add_argument("ppi_output", help="PPi file to be integrated in neXtProt")
parser.add_argument("interaction_mapping_output", help="interaction mapping file to be integrated in neXtProt")
parser.add_argument("step", help="[psimi_merged, psimi_unmerged]: 'psimi_merged' to concatenate psimi, 'psimi_unmerged' new psimi new line", choices=["psimi_merged", "psimi_unmerged"])
args = parser.parse_args()

logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', filename='main.log', level=logging.DEBUG)

###The aim of this script is to create one output file with PPi to be integrated and one output file with interaction mapping information to be integrated

## INPUTS: ##
##-The third-version file (3/4) with PPi and propagated PPi thanks to identical protein sequence
##-[psimi_merged, psimi_unmerged] to concatenate psimi when other information are identical

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
    output_line += "\tENYO" ### tout ce qui nest pas ppi["interactor1_accession"] on peut supprimer
    output_line += "\tcurated"
    output_line += "\tPROTEIN"
    output_line += "\tNX_"+ppi["interactor1_accession"]
    output_line += "\tECO:0000353"
    output_line += "\tGOLD"
    output_line += "\tPubMed"
    output_line += "\tpublication"
    output_line += "\tENYO"
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

    ##for VH-PPI modify here
    # if (ppi["ppi_type"] == "vh"):
    #     output_line += "\t"+ppi["Uniprot_viral_species"]+"\t"+ppi["FTId"]+"\t"+ppi["Chain_name"]
    #     output_line += "\t"+ppi["interAC_propagation"]
    # elif (ppi["ppi_type"] == "hh"):
    #     output_line += "\t\t\t"
    #     output_line += "\t"+ppi["interAC_propagation"]

    intmap_list.append(output_line)


##Function to store every interaction mapping information
def integration_formatting():
    logging.info("Starts integration formatting")
    ppi_result = open(args.ppi_output,"w")
    intmap_result = open(args.interaction_mapping_output, "w")

    ##write the header in PPi output file
    ppi_header = "reference_accession\tpsimi_id"
    ppi_header += "\tnextprot_accession\tgene_name" #interactant1
    ppi_header += "\tbiological_object_accession\tbiological_object_name" #interactant2
    ##for VH-PPI modify here
    #ppi_header += "\tUniprot_viral_species\tFTId\tChain_name"
    #ppi_header += "\tinterAC_propagation"

    #extra values
    ppi_header += "\tannotation_category"
    ppi_header += "\tassigned_by"
    ppi_header += "\tassigment_method"
    ppi_header += "\tbiological_object_type"
    ppi_header += "\tentry_accession"
    ppi_header += "\tevidence_code"
    ppi_header += "\tevidence_quality"
    ppi_header += "\treference_database"
    ppi_header += "\tresource_type"
    ppi_header += "\tsource"
    ppi_header += "\tannot_source_accession"

    ppi_result.write(ppi_header.upper()+"\n")

    ##write the header in interaction mapping output file
    intmap_header = "pmid\tpsimi_id"
    intmap_header += "\tinteractor1_accession\tinteractor1_name"
    intmap_header += "\tinteractor2_accession\tinteractor2_name"
    intmap_header += "\tinteractor_mapping_sequence\tinteractor_isoform_accession"
    intmap_header += "\tinteractor_occurrence_start\tinteractor_occurrence_stop\tinteractor_occurrence_identity"
    ##for VH-PPI modify here
    # intmap_header += "\tUniprot_viral_species\tFTId\tChain_name"
    # intmap_header += "\tinterAC_propagation"
    intmap_result.write(intmap_header.upper()+"\n")
    
    ppi_list=[]
    intmap_list=[]

    with open(args.propagated_ppi_file) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        ##fill lists of ppi and interaction mapping
        for ppi in reader:
            store_ppi(ppi, ppi_list)
            if ( (ppi["interactor1_isoform_accession"] != "") or (ppi["interactor2_isoform_accession"] != "") ):
                store_interaction_mapping(ppi, intmap_list)

        ##write PPi
        ppi_list = list(set(ppi_list)) ##save one copy of the duplicated PPi, and sort them
        ppi_list.sort(key=lambda x: x.split("\t")[4]) #sort by int2
        ppi_list.sort(key=lambda x: x.split("\t")[2]) #sort by int1
        ppi_list.sort(key=lambda x: x.split("\t")[0]) #sort by pmid, mandatory to compare lines by lines
        #print(ppi_list)

        cpt = 0

        if (args.step == "psimi_unmerged"):
            for uniq_ppi in ppi_list: ##write the full list
                cpt += 1
                ppi_result.write(uniq_ppi+"\t"+str(cpt)+"\n")


        if (args.step == "psimi_merged"):

            prevppi = ppi_list[0] ##first element from the list
            prevkey_list = prevppi.split("\t")
            
            prevpmid = prevkey_list[0]
            prevpsimi = prevkey_list[1]
            prevint1 = prevkey_list[2]
            prevname1 = prevkey_list[3]
            prevint2 = prevkey_list[4]
            prevname2 = prevkey_list[5]
            previsoform = prevkey_list[10]

            prevkey = prevpmid+"_"+prevint1+"_"+prevint2
            multiple_psimi = prevpsimi
            #print(prevkey)

            for uniq_ppi in ppi_list[1:]: ##every other element than the first element from the list
                key_list = uniq_ppi.split("\t")
                mypmid = key_list[0]
                mypsimi = key_list[1]
                myint1 = key_list[2]
                myname1 = key_list[3]
                myint2 = key_list[4]
                myname2 = key_list[5]
                myisoform = key_list[10]

                mykey = mypmid+"_"+myint1+"_"+myint2
                #print(mykey)

                if ( prevkey != mykey):
                    cpt += 1
                    myppi = prevpmid+"\t"+multiple_psimi+"\t"+prevint1+"\t"+prevname1+"\t"+prevint2+"\t"+prevname2
                    myppi += "\tbinary-interaction\tENYO\tcurated\tPROTEIN\t"+myisoform
                    myppi += "\tECO:0000353\tGOLD\tPubMed\tpublication\tENYO"
                    ppi_result.write(myppi+"\t"+str(cpt)+"\n")
                    
                    prevpmid = mypmid
                    prevpsimi = mypsimi
                    prevint1 = myint1
                    prevname1 =myname1
                    prevint2 = myint2
                    prevname2 = myname2
                    previsoform = myisoform

                    prevkey = mykey
                    multiple_psimi = mypsimi
                    
                else:
                    multiple_psimi += ", "+mypsimi #one ppi one pmid but several psimi

                ##if last element of the list write it anyway
                if uniq_ppi == ppi_list[-1]:
                    cpt += 1
                    myppi = mypmid+"\t"+multiple_psimi+"\t"+myint1+"\t"+myname1
                    myppi += "\t"+myint2+"\t"+myname2
                    myppi += "\tbinary-interaction\tENYO\tcurated\tPROTEIN\t"+myisoform
                    myppi += "\tECO:0000353\tGOLD\tPubMed\tpublication\tENYO"
                    ppi_result.write(myppi+"\t"+str(cpt)+"\n")

        ##write intmap
        intmap_list = list(set(intmap_list)) ##save only one copy of homodimer self-interaction (interactor1 and interactor2 have both the same intmap region)
        intmap_list.sort(key=lambda x: x.split("\t")[7])
        intmap_list.sort(key=lambda x: x.split("\t")[6])
        intmap_list.sort(key=lambda x: x.split("\t")[4])
        intmap_list.sort(key=lambda x: x.split("\t")[2])
        intmap_list.sort(key=lambda x: x.split("\t")[0])

        if (args.step == "psimi_unmerged"):
            for uniq_intmap in intmap_list:
                intmap_result.write(uniq_intmap+"\n")

        if (args.step == "psimi_merged"):
            previntmap = intmap_list[0] #first_line
            prevkey_list = previntmap.split("\t")

            prevpmid = prevkey_list[0]
            prevpsimi = prevkey_list[1]
            prevint1 = prevkey_list[2]
            prevname1 = prevkey_list[3]
            prevint2 = prevkey_list[4]
            prevname2 = prevkey_list[5]
            prevmapseq = prevkey_list[6]
            previsoform = prevkey_list[7]
            prevoccstart = prevkey_list[8]
            prevoccstop = prevkey_list[9]
            prevoccidentity = prevkey_list[10]
        
            prevkey = prevpmid+"_"+prevint1+"_"+prevint2+"_"+previsoform+"_"+prevmapseq
            multiple_psimi = prevpsimi
            #print(prevkey)

            for uniq_intmap in intmap_list[1:]:
                key_list = uniq_intmap.split("\t")
                mypmid = key_list[0]
                mypsimi = key_list[1]
                myint1 = key_list[2]
                myname1 = key_list[3]
                myint2 = key_list[4]
                myname2 = key_list[5]
                mymapseq = key_list[6]
                myisoform = key_list[7]
                myoccstart = key_list[8]
                myoccstop = key_list[9]
                myoccidentity = key_list[10]
                mykey = mypmid+"_"+myint1+"_"+myint2+"_"+myisoform+"_"+mymapseq
                #print(mykey)

                if ( prevkey != mykey):
                    myintmap = prevpmid+"\t"+multiple_psimi+"\t"+prevint1+"\t"+prevname1
                    myintmap += "\t"+prevint2+"\t"+prevname2+"\t"+prevmapseq+"\t"+previsoform
                    myintmap += "\t"+prevoccstart+"\t"+prevoccstop+"\t"+prevoccidentity
                    intmap_result.write(myintmap+"\n")

                    prevpmid = mypmid
                    prevpsimi = mypsimi
                    prevint1 = myint1
                    prevname1 = myname1
                    prevint2 = myint2
                    prevname2 = myname2
                    prevmapseq = mymapseq
                    previsoform = myisoform
                    prevoccstart = myoccstart
                    prevoccstop = myoccstop
                    prevoccidentity = myoccidentity
                    
                    prevkey = mykey
                    multiple_psimi = mypsimi

                else:
                    multiple_psimi += ", "+mypsimi

                if uniq_intmap == intmap_list[-1]:
                    myintmap = mypmid+"\t"+multiple_psimi+"\t"+myint1+"\t"+myname1
                    myintmap += "\t"+myint2+"\t"+myname2+"\t"+mymapseq+"\t"+myisoform
                    myintmap += "\t"+myoccstart+"\t"+myoccstop+"\t"+myoccidentity
                    intmap_result.write(myintmap+"\n")
                            
    ppi_result.close()
    intmap_result.close()
    logging.info("Integration formatting ends")

####################
# MAIN #############
####################

if __name__ == "__main__":
    integration_formatting() #sys.argv[1], sys.argv[2], sys.argv[3]
#./integration_formatting.py integre_vinland-hh-2020-01-22.tsv_one_iso_propagated_symmetric myppi2 myintmap2 psimi_merged
#./integration_formatting.py integre_vinland-hh-2020-01-22.tsv_one_iso_propagated_symmetric myppi1 myintmap1 psimi_unmerged
