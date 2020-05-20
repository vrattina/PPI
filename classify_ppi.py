#AUTHOR: RATTINA Vimel - 2019/11/28 - SIB & Enyo Pharma

import sys #I/O files
import os, errno #create folder
import csv #parse tsv
import logging #log file
import pickle #store and use python object

logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', filename='main.log', level=logging.DEBUG)

###The aim of this script is to create a file to be integrated in neXtProt with no gene name and/or sequence changes and a file to be curated due to gene name and/or sequence changes

## INPUTS: ##
##-The flatten PPi file from ENYO (output of flatten_enyo_json.py),
##-The dictionnary object file with viral protein information about gene name and sequence from 2 versions of uniprot and if there is changes or not.
##-The dictionnary object file with human protein information about gene name and sequence from 2 versions of uniprot and if there is changes or not.
##-An output file name for the file to be integrated in neXtProt
##-An output file name for the file to be recurated by enyo & neXtProt

## OUTPUTS: ##
##-A first-version file (1/4) to be integrated by neXtProt with no gene name or sequence changes between 2 versions of uniprot i.e. the flatten PPi with in addition viral species, chain name and FTId
##-A file to be curated by enyo & nextprot with gene name and/or sequence changes between 2 versions of uniprot

###################
# FUNCTIONS #######
###################

###Function writting the output file ready to be integrated
def no_changes_for_sib(ppi, interactor1_ID, interactor2_ID, HP_dict, HorV_P_dict, outputfile):
    ##if changes status is available, the AC get downloaded
    if ( ("Genename_changes" in HP_dict[interactor1_ID]) and ("Genename_changes" in HorV_P_dict[interactor2_ID])):
        ##if no gene name or sequence changed
        if ( (HorV_P_dict[interactor2_ID]["Genename_changes"] == "no") and (HorV_P_dict[interactor2_ID]["Sequence_changes"] == "no") and (HP_dict[interactor1_ID]["Genename_changes"] == "no") and (HP_dict[interactor1_ID]["Sequence_changes"] == "no")):
            common_core = ppi["ppi_type"]+"\t"+ppi["stable_id"]+"\t"+ppi["pmid"]+"\t"+ppi["psimi_id"]
            common_core += "\t"+ppi["interactor1_accession"]+"\t"+ppi["interactor1_name"]+"\t"+ppi["interactor1_start"]+"\t"+ppi["interactor1_stop"]
            common_core += "\t"+ppi["interactor2_accession"]+"\t"+ppi["interactor2_name"]+"\t"+ppi["interactor2_start"]+"\t"+ppi["interactor2_stop"]
            common_core += "\t"+str(ppi["interactor1_mapping_sequence"])+"\t"+str(ppi["interactor1_isoform_accession"])+"\t"+str(ppi["interactor1_occurrence_start"])+"\t"+str(ppi["interactor1_occurrence_stop"])+"\t"+str(ppi["interactor1_occurrence_identity"])
            common_core += "\t"+str(ppi["interactor2_mapping_sequence"])+"\t"+str(ppi["interactor2_isoform_accession"])+"\t"+str(ppi["interactor2_occurrence_start"])+"\t"+str(ppi["interactor2_occurrence_stop"])+"\t"+str(ppi["interactor2_occurrence_identity"])

            ##if VH PPi type every description will be associated to a viral species
            if (ppi["ppi_type"] == "vh"):
                common_core += "\t"+HorV_P_dict[interactor2_ID]["Uniprot_viral_species"]

                ##when chain name and FTId available write it otherwise write None
                if HorV_P_dict[interactor2_ID].has_key("Chain_name"):
                    common_core += "\t"+HorV_P_dict[interactor2_ID]["FTId"]+"\t"+HorV_P_dict[interactor2_ID]["Chain_name"]
                else:
                    common_core += "\tNone\tNone"

            elif (ppi["ppi_type"] == "hh"):
                common_core += "\t\t\t"

            outputfile.write(common_core+"\n")


###Function writting the output file to be curated because of changes in gene name or sequence between the 2 UniProt versions
def changes_to_cure(ppi, interactor1_ID, interactor2_ID, HP_dict, HorV_P_dict, outputfile):
    ##if changes status is available, the AC get downloaded
    if ( ("Genename_changes" in HP_dict[interactor1_ID]) and ("Genename_changes" in HorV_P_dict[interactor2_ID])):

        ##if changes repertoried in HP, print the pmid URL, gene name and sequence from the two UniProt versions
        if ( (HP_dict[interactor1_ID]["Genename_changes"] == "yes") or (HP_dict[interactor1_ID]["Sequence_changes"] == "yes")):
            print(interactor1_ID+"\tchanges")
            debug = ppi["ppi_type"]+"\t"+ppi["stable_id"]+"\thttps://www.ncbi.nlm.nih.gov/pubmed/"+ppi["pmid"]
            debug += "\t"+ppi["interactor1_accession"]+"\t"+ppi["interactor1_name"]
            debug += "\t"+HP_dict[interactor1_ID]["Genename_changes"]+"\t"+HP_dict[interactor1_ID]["Sequence_changes"]
            debug += "\t"+HP_dict[interactor1_ID]["Genename_ENYO"]+"\t"+HP_dict[interactor1_ID]["Genename_neXtProt"]
            debug += "\t"+HP_dict[interactor1_ID]["Sequence_ENYO"]+"\t"+HP_dict[interactor1_ID]["Sequence_neXtProt"]
            debug += "\t"+str(ppi["interactor1_isoform_accession"])+"\tHomo sapiens" ##Here the species is Homo sapiens

            outputfile.write(debug+"\n")

        ##If changes repertoried in interactor2 so VP or HP depending on the PPi type
        if ( (HorV_P_dict[interactor2_ID]["Genename_changes"] == "yes") or (HorV_P_dict[interactor2_ID]["Sequence_changes"] == "yes")):
            print(interactor2_ID+"\tchanges")
            debug = ppi["ppi_type"]+"\t"+ppi["stable_id"]+"\thttps://www.ncbi.nlm.nih.gov/pubmed/"+ppi["pmid"]
            debug += "\t"+ppi["interactor2_accession"]+"\t"+ppi["interactor2_name"]
            debug += "\t"+HorV_P_dict[interactor2_ID]["Genename_changes"]+"\t"+HorV_P_dict[interactor2_ID]["Sequence_changes"]
            debug += "\t"+HorV_P_dict[interactor2_ID]["Genename_ENYO"]+"\t"+HorV_P_dict[interactor2_ID]["Genename_neXtProt"]
            debug += "\t"+HorV_P_dict[interactor2_ID]["Sequence_ENYO"]+"\t"+HorV_P_dict[interactor2_ID]["Sequence_neXtProt"]
            debug += "\t"+str(ppi["interactor2_isoform_accession"])

            if (ppi["ppi_type"] == "vh"):
                debug += "\t"+HorV_P_dict[interactor2_ID]["Uniprot_viral_species"]
            elif (ppi["ppi_type"] == "hh"):
                debug += "\tHomo sapiens"

            outputfile.write(debug+"\n")


###Function filling the output file to be curated because no information was found
def no_information(ppi, interactor1_ID, interactor2_ID, HP_dict, HorV_P_dict, outputfile):
    ##if Genename_changes key does not exist, a problem occurs during download step
    ##if no changes in HP AC
    if ("Genename_changes" not in HP_dict[interactor1_ID]):
        print(interactor1_ID+"\tdoes not exist")
        debug = ppi["ppi_type"]+"\t"+ppi["stable_id"]+"\thttps://www.ncbi.nlm.nih.gov/pubmed/"+ppi["pmid"]
        debug += "\t"+ppi["interactor1_accession"]+"\t"+ppi["interactor1_name"]
        debug += "\tyes\tyes\tNone\tNone\tNone\tNone" ##no gene name and sequence so put yes for both changes and None because no data were retrieved
        debug += "\t"+str(ppi["interactor1_isoform_accession"]+"\thuman_failed") ##add a species name such as failed to highlight that the data could not be retrieved
        outputfile.write(debug+"\n")

    ##If no changes in interactor2 so VP AC or HP AC
    if ("Genename_changes" not in HorV_P_dict[interactor2_ID]):
        print(interactor2_ID+"\tdoes not exist")
        debug = ppi["ppi_type"]+"\t"+ppi["stable_id"]+"\thttps://www.ncbi.nlm.nih.gov/pubmed/"+ppi["pmid"]
        debug += "\t"+ppi["interactor2_accession"]+"\t"+ppi["interactor2_name"]
        debug += "\tyes\tyes\tNone\tNone\tNone\tNone"

        if (ppi["ppi_type"] == "vh"):
            debug += "\tNone\tvirus_failed"
        elif (ppi["ppi_type"] == "hh"):
            debug += "\tNone\thuman_failed"
        outputfile.write(debug+"\n")


##Function creating the output files
def classify_ppi(flatfile, VP_dict_file, HP_dict_file, outputfile_raw_integrate, outputfile_recurate):
    logging.info("Starts classifying PPi")
    output_sib = open(outputfile_raw_integrate,"w") 
    output_cure = open(outputfile_recurate,"w")

    with open(VP_dict_file, "rb") as viral_file:
        VP_dict = pickle.loads(viral_file.read())

    with open(HP_dict_file, "rb") as human_file:
        HP_dict = pickle.loads(human_file.read())
    
    ##write the header in both output files
    header_integrate = "ppi_type\tstable_id\tpmid\tpsimi_id"

    header_integrate += "\tinteractor1_accession\tinteractor1_name\tinteractor1_start\tinteractor1_stop"
    header_integrate += "\tinteractor2_accession\tinteractor2_name\tinteractor2_start\tinteractor2_stop"

    header_integrate += "\tinteractor1_mapping_sequence\tinteractor1_isoform_accession"
    header_integrate += "\tinteractor1_occurrence_start\tinteractor1_occurrence_stop\tinteractor1_occurrence_identity" 

    header_integrate += "\tinteractor2_mapping_sequence\tinteractor2_isoform_accession"
    header_integrate += "\tinteractor2_occurrence_start\tinteractor2_occurrence_stop\tinteractor2_occurrence_identity" 

    header_integrate += "\tUniprot_viral_species\tFTId\tChain_name\n"

    output_sib.write(header_integrate)

    header_cure = "ppi_type\tstable_id\tpmid"
    header_cure += "\tinteractor_accession\tinteractor_protein_name\tGenename_changes\tSequence_changes"
    header_cure += "\tGenename_ENYO\tGenename_neXtProt\tSequence_ENYO\tSequence_neXtProt"
    header_cure += "\tinteractor_isoform_AC\tspecies\n"

    output_cure.write(header_cure)

    with open(flatfile) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for ppi in reader:
            
            if (ppi["ppi_type"] == "vh"):
                ##recreate the keys, polyprotein name for viral protein and with isoform for human protein
                interactor1_ID = ""
                if (ppi["interactor1_mapping_sequence"] != ""): #None
                    interactor1_ID = ppi["interactor1_isoform_accession"]
                else:
                    interactor1_ID = ppi["interactor1_accession"]

                interactor2_ID = ppi["interactor2_accession"]+"_"+ppi["interactor2_name"]+"_"+ppi["interactor2_start"]+"_"+ppi["interactor2_stop"]
                ##function to write in the output files
                no_changes_for_sib(ppi, interactor1_ID, interactor2_ID, HP_dict, VP_dict, output_sib)
                changes_to_cure(ppi, interactor1_ID, interactor2_ID, HP_dict, VP_dict, output_cure)
                no_information(ppi, interactor1_ID, interactor2_ID, HP_dict, VP_dict, output_cure)


            if (ppi["ppi_type"] == "hh"):
                interactor1_ID = ""
                if (ppi["interactor1_mapping_sequence"] != ""): #None
                    interactor1_ID = ppi["interactor1_isoform_accession"]
                else:
                    interactor1_ID = ppi["interactor1_accession"]
                    #print VP_dict[interactor1_ID]["Vprot_accession"]
                
                interactor2_ID = ""
                if (ppi["interactor2_mapping_sequence"] != ""): #None
                    interactor2_ID = ppi["interactor2_isoform_accession"]
                else:
                    interactor2_ID = ppi["interactor2_accession"]
                    #print VP_dict[interactor1_ID]["Vprot_accession"]

                no_changes_for_sib(ppi, interactor1_ID, interactor2_ID, HP_dict, HP_dict, output_sib)
                changes_to_cure(ppi, interactor1_ID, interactor2_ID, HP_dict, HP_dict, output_cure)
                no_information(ppi, interactor1_ID, interactor2_ID, HP_dict, HP_dict, output_cure)

    output_cure.close()
    output_sib.close()
    logging.info("PPi classifying ends")

####################
# MAIN #############
####################

if len(sys.argv) < 6:
    sys.exit('Usage: %s <flat_file> <VP_dict> <HP_dict> <outputfile_raw_integrate> <outputfile_recurate>\n<flat_file>: flatten json file from ENYO Pharma\n<VP_dict>: dictionnary of viral protein from update_ac_uniprot_information function\n<HP_dict>: dictionnary of human protein from update_ac_uniprot_information function\n<outputfile_raw_integrate>: first version (/4) of the file to be integrated in neXtProt (no gene name or sequence changes) with the same columns than flat_file\n<outputfile_recurate>: file to be recurate because of gene name and/or sequence changes or AC missing' % sys.argv[0])

if __name__ == "__main__":
    classify_ppi(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
#./classify_ppi.py vinland_flatfile.tsv results/2019_01/ results/2019_09/ theoutput_sib.tsv theoutput_recurate.tsv &> theoutput2.log
