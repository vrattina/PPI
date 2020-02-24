#!/usr/bin/python
#AUTHOR: RATTINA Vimel - 2020/01/08 - SIB & Enyo Pharma

import sys #I/O files
import os, errno #create folder
import csv #parse tsv
import logging #log file
import subprocess #run perl script on python
import re #regular expression
import pickle #store and use python object

logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', filename='main.log', level=logging.DEBUG)

### This script has DEPENDENCIES: varsplic.pl and swiss knife (http://swissknife.sourceforge.net/)

##The aim of this script is to create two distinct dictionnaries of AC (one for human proteins and one for viral proteins) with their gene name(s), sequence, if there is gene name or sequence changes and additionnal information for viral proteins such as the viral species, the corresponding chain name and feature id (FTId).

## INPUTS: ##
##-The flatten PPi file from ENYO (output of flatten_enyo_json.py),
##-Pathway if the result folder which contains the two UniProt version repositories,
##-the ENYO UniProt version,
##-And the neXtProt UniProt version.

## OUTPUTS: ##
##-A human protein dictionnary object with gene name, sequence, gene name changes and sequence changes,
##-A viral protein dictionnary object with protein name, sequence, genename changes, sequence changes, viral species, the corresponding chain name and feature id FTId.

###################
# FUNCTIONS #######
###################

###Function to fill the human protein dictionnary with the corresponding information from UniProt
def hp_filler(ppi, HP_dict, interactor):
    ##variable initialization with the wanted interactor name
    mapping_sequence=interactor+"_mapping_sequence"
    accession=interactor+"_accession"
    name=interactor+"_name"
    start=interactor+"_start"
    stop=interactor+"_stop"
    isoform_accession=interactor+"_isoform_accession"
    occurrence_start=interactor+"_occurrence_start"
    occurrence_stop=interactor+"_occurrence_stop"
    occurrence_identity=interactor+"_occurrence_identity"

    ##Fill human protein HP dictionnary
    if (ppi[mapping_sequence] == "") :
        HP_dict[ppi[accession]]={}
        HP_dict[ppi[accession]].update(Hprot_accession=ppi[accession], Hgene_name=ppi[name], Hprot_start=ppi[start], Hprot_stop=ppi[stop])
    else :
        ##If interaction mapping information is available, use the correct isoform such as key instead of the canonical isoform AC by default
        HP_dict[ppi[isoform_accession]]={}
        HP_dict[ppi[isoform_accession]].update(Hprot_accession=ppi[isoform_accession], Hgene_name=ppi[name], Hprot_start=ppi[start], Hprot_stop=ppi[stop], mapping_sequence=ppi[mapping_sequence], occurrence_start=ppi[occurrence_start], occurrence_stop=ppi[occurrence_stop], occurrence_identity=ppi[occurrence_identity])
        #print ppi['isoform_accession']


###Function creating two dictionnaries, one storing viral proteins (if VH-PPi available) and one storing human proteins
def vp_hp_in_dictionnary(flatfile):
    VP, HP = {}, {}
    with open(flatfile) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')

        for ppi in reader:

            if (ppi['ppi_type'] == "vh"):
                ##interactor1 is the human protein, use hp_filler to fill HP dictionnary
                hp_filler(ppi, HP, "interactor1")

                ##creation of a new viral protein entry such as dictionnary key to treat polyprotein (same AC however several proteins)
                VP_mature_name=ppi['interactor2_accession']+"_"+ppi['interactor2_name']+"_"+ppi['interactor2_start']+"_"+ppi['interactor2_stop'] 
                ##fill viral protein (VP) dictionnary
                VP[VP_mature_name]={}
                VP[VP_mature_name].update(Vprot_accession=ppi['interactor2_accession'], Vprot_name=ppi['interactor2_name'], Vprot_start=ppi['interactor2_start'], Vprot_stop=ppi['interactor2_stop'])

                
            elif (ppi['ppi_type'] == "hh"):
                ##in HH PPi type, both interactors are human proteins so use hp_filler twice
                hp_filler(ppi, HP, "interactor1")
                hp_filler(ppi, HP, "interactor2")

    return VP, HP


###Function treating the dictionnaries to retrieve every AC by removing the duplicated data
def uniq_ac(VP_dict, HP_dict):
    VP_AC=[]
    HP_AC=[]
    ##For VP do not take the polyprotein key but the AC value
    for AC in VP_dict.values():
        VP_AC.append(AC["Vprot_accession"])
    VP_AC = list(set(VP_AC)) #uniq AC

    ##For HP, retrive keys because they own AC and isoform
    HP_AC += HP_dict.keys() #merge lists
    HP_AC = list(set(HP_AC))

    return VP_AC, HP_AC


###Function retrieving gene names and protein sequence from a AC.txt file (only canonical isoform is available)
def get_genename_sequence_canonical(ACpath):
    gene_name=[]
    next_line_is_sequence=False
    sequence=""
    txt = open(ACpath, "r")
    for i in txt.readlines():
        i = i.rstrip()

        if i.startswith("GN   "):
            gene_name += i.replace('ORF ','ORF').replace('GN','').replace('Name=','').replace('Synonyms=','').replace(',','').replace('OrderedLocusNames=','').replace('ORFNames=','').replace(';','').replace('-','').replace('_','').replace('/','|').upper().split() ##If ORF 63, transform to ORF63

        if (next_line_is_sequence == True): ##Up to the last line // for end of file
            sequence+=i

        if i.startswith("SQ   SEQUENCE"):
            next_line_is_sequence = True

    ##Remove digit name
    ##Remove name containing : e.g. ECO: etc
    intermediar_genename = [ x for x in gene_name if ( (":" not in x) and (not x.isdigit()) ) ] 
    gene_name = set(intermediar_genename) ##Type is set to permit intersection of gene names between two UniProt versions
    fasta = sequence.replace(" ", "").replace("//", "")

    return gene_name, fasta


###Function retrieving protein sequence for one alternate isoform chosen
###varsplic.pl and swiss knife dependencies
def get_sequence_isoform(ACpath, isoformAC, YYYY_MM):
    varsplic_script=""

    ##according to UniProt version use the corresponding varsplic.pl script using the adapted swissknife version
    if (YYYY_MM >= "2019_11"):
        varsplic_script = "/scratch/Softwares/swissknife_1.79/lib/varsplic.pl"
    elif (YYYY_MM < "2019_11"):
        varsplic_script = "/scratch/Softwares/swissknife_1.78/lib/varsplic.pl"

    job = subprocess.Popen(["perl", varsplic_script, "-input", ACpath, "-fasta"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    alt_iso_fasta, errors = job.communicate()

    ### Regexp method
    fasta_pattern=r"(>)sp\|(" + re.escape(isoformAC) + r")\|[^\n]+([^>]+)"
    gene_name_sequence = re.search(fasta_pattern, alt_iso_fasta, re.DOTALL)
    
    if gene_name_sequence is None:
        logging.warning(isoformAC+" regular expression does not work, it is probable that the alternate isoform does not exist please check manually: perl /path/swissknife_version/varsplic.pl -input "+ACpath+" -fasta")
        fasta = "None"
    else:
        fasta = gene_name_sequence.group(3).replace("\n", "")
        #fasta = gene_name_sequence.group(1)+gene_name_sequence.group(2)+"\n"+gene_name_sequence.group(3).replace("\n", "")
    #print isoformAC+"\t"+ACpath+"\t"+alt_iso_fasta
    return fasta


###Function adding information in VP or HP dictionnary concerning gene names and protein sequence from both UniProt version and if there is changes between versions or not
def update_changes_status(P_dict, isoform, sequence_changes, genename_changes, name_enyo, name_nextprot, sequence_enyo, sequence_nextprot):
    if ("Vprot_accession" in P_dict.values()[0]):
        for i in P_dict.values():
            if i["Vprot_accession"] == isoform:
                i["Sequence_changes"] = sequence_changes
                i["Genename_changes"] = genename_changes
                i["Genename_ENYO"] = name_enyo
                i["Genename_neXtProt"] = name_nextprot
                i["Sequence_ENYO"] = sequence_enyo
                i["Sequence_neXtProt"] = sequence_nextprot

    if ("Hprot_accession" in P_dict.values()[0]):
        for i in P_dict.values():
            if i["Hprot_accession"] == isoform:
                i["Sequence_changes"] = sequence_changes
                i["Genename_changes"] = genename_changes
                i["Genename_ENYO"] = name_enyo
                i["Genename_neXtProt"] = name_nextprot
                i["Sequence_ENYO"] = sequence_enyo
                i["Sequence_neXtProt"] = sequence_nextprot


###Function comparing gene names and sequence of one AC from two UniProt versions and updating the protein dictionnary
def compare_ac_versions(isoform_name, P_dict, outputfolder, YYYY_MM_enyo, YYYY_MM_nextprot):
    if '-' in isoform_name:
        canonical_name=isoform_name.split('-')[0]
    else:
        canonical_name = isoform_name
    #print isoform_name+"\t"+canonical_name

    uniprot_enyo_folder = outputfolder+"/"+YYYY_MM_enyo
    uniprot_nextprot_folder = outputfolder+"/"+YYYY_MM_nextprot
    path_enyo = uniprot_enyo_folder+"/"+canonical_name+".txt"
    path_nextprot = uniprot_nextprot_folder+"/"+canonical_name+".txt"
    
    if os.path.exists(path_enyo) and os.path.exists(path_nextprot):
        ##retrieve sequence and gene names from two UniProt version
        name_enyo, sequence_enyo = get_genename_sequence_canonical(path_enyo)
        name_nextprot, sequence_nextprot = get_genename_sequence_canonical(path_nextprot)

        clean_name_enyo = ','.join(list(name_enyo))
        clean_name_nextprot = ','.join(list(name_nextprot))
        
        ##if alternative isoform (e.g. AC-1), retrieve the correct sequence
        if '-' in isoform_name:
            sequence_enyo = get_sequence_isoform(path_enyo, isoform_name, YYYY_MM_enyo)
            sequence_nextprot = get_sequence_isoform(path_nextprot, isoform_name, YYYY_MM_nextprot)
            
        ##check if there is changes or not between the 2 UniProt version
        sequence_changes = ""
        genename_changes = ""
        if (sequence_enyo != sequence_nextprot):
            logging.warning(isoform_name+" protein sequence changes between "+uniprot_enyo_folder+" and "+uniprot_nextprot_folder)
            sequence_changes="yes"
        else:
            sequence_changes="no"

        ##if not even one synonym or one gene name shared, and both names are not identical (e.g. both empty)
        if (not (name_enyo.intersection(name_nextprot))) and (name_enyo != name_nextprot) :
            logging.warning(canonical_name+" gene name changes between "+uniprot_enyo_folder+" and "+uniprot_nextprot_folder)
            genename_changes="yes"
        else:
            genename_changes="no"

        ##add information about gene name, sequence and the changes status in the corresponding dictionnary
        update_changes_status(P_dict, isoform_name, sequence_changes, genename_changes, clean_name_enyo, clean_name_nextprot, sequence_enyo, sequence_nextprot)
    #print P_dict


###Function adding viral species information and chain name with feature id (FTId) if corresponding positions
def get_species_chain_ftid(ACpath, VP_dict):
    txt = open(ACpath, "r")
    ftid_line=False ##FTId is in few lines after chain name
    species = "" ##Initiate because some protein can have several OS description
    for i in txt.readlines():
        i = i.rstrip()
        
        if i.startswith("OS   "):
            species += i.split("OS   ")[1].split(".")[0]
            VP_dict['Uniprot_viral_species'] = species

        if (ftid_line == True):
            ##double if because the FTId is in few lines after chain name (not always the next one)

            ##if UniProt version < 2019_11, FTId and chain name in the same line, then break
            if "FTId=" in i: 
                ftid = i.split("=")[1].split('.')[0]
                VP_dict['FTId']=ftid
                VP_dict['Chain_name']=chain_name
                #print VP_dict['Vprot_name']+"\t"+chain_name+" "+ftid
                break

            ##if UniProt version >= 2019_11, FTId and chain name in different lines, the last information is the FTId
            if "note=" in i: #1..50
                chain_name = i.split('=')[1].replace('"','')
                VP_dict['Chain_name']=chain_name
            if "id=" in i:
                ftid = i.split("=")[1].split('.')[0]
                VP_dict['FTId']=ftid
                break
            
        if i.startswith("FT   CHAIN"):
            ftchain = re.sub(r"  +", "\t", i).split('\t')
            chain_start=""
            chain_stop=""
            if (len(ftchain) >= 5): #if position space separated e.g. 1 50, UniProt version < 2019_11
                chain_start = ftchain[2]
                chain_stop = ftchain[3]
                chain_name = ftchain[4].split('.')[0]

            elif ".." in ftchain[2]: #if position .. separated e.g. 1..50, UniProt version >= 2019_11
                chain_start=ftchain[2].split('..')[0]
                chain_stop=ftchain[2].split('..')[1]

            if ((VP_dict['Vprot_start'] == chain_start) and (VP_dict['Vprot_stop'] == chain_stop)): 
                
                ftid_line = True


###Function returning two dictionnaries with gene name, sequence, if any changes and additionnal information for viral proteins (viral species and chain name with FTId if available)
def update_ac_uniprot_information(flatfile, outputfolder, YYYY_MM_enyo, YYYY_MM_nextprot, vp_dict_file, hp_dict_file):
    logging.info("Starts AC comparison between "+YYYY_MM_enyo+ " and "+YYYY_MM_nextprot+" to check if any gene name or sequence changes")

    ##Put PPi information in a dictionnary
    VP_dict, HP_dict = vp_hp_in_dictionnary(flatfile)
    ##Retrieve two lists of AC for human and viral proteins
    VP_AC, HP_AC = uniq_ac(VP_dict, HP_dict)

    ##Compare two UniProt versions
    logging.info("Comparing AC between "+YYYY_MM_enyo+ " and "+YYYY_MM_nextprot+" for viral proteins")
    for i in VP_AC:
        compare_ac_versions(i, VP_dict, outputfolder, YYYY_MM_enyo, YYYY_MM_nextprot)
    logging.info("Viral protein comparison ends")

    logging.info("Comparing AC between "+YYYY_MM_enyo+ " and "+YYYY_MM_nextprot+" for human proteins")
    for i in HP_AC:
        compare_ac_versions(i, HP_dict, outputfolder, YYYY_MM_enyo, YYYY_MM_nextprot)
    logging.info("Human protein comparison ends")

    ##Add viral species + chain name and FTId if available for VP
    logging.info("Adding species name for viral proteins and chain name and PRO id if corresponding to the positions")
    for i in VP_dict.values():
        #print i["Vprot_accession"]+" "+i['Vprot_start']+" "+i['Vprot_stop']
        AC_nextprot = outputfolder+"/"+YYYY_MM_nextprot+"/"+i["Vprot_accession"]+".txt" ##Version of neXtProt, data will be integrated in it
        if os.path.exists(AC_nextprot):
            #print AC_nextprot
            get_species_chain_ftid(AC_nextprot, i)
    logging.info("Viral protein information completion ends")

    with open(vp_dict_file, "w") as viral_file:
        pickle.dump(VP_dict, viral_file)

    with open(hp_dict_file, "w") as human_file:
        pickle.dump(HP_dict, human_file)

    logging.info("AC comparison between "+YYYY_MM_enyo+ " and "+YYYY_MM_nextprot+" ends")
    
####################
# MAIN #############
####################

if len(sys.argv) < 7:
    sys.exit('Usage: %s <flat_file> <output_folder> <YYYY_MM_enyo> <uniprotv_nextprot_folder> <vp_dict_file> <hp_dict_file>\n<flat_file>: flatten json file from ENYO Pharma\n<outputfolder>: pathway of the result folder which contains the two repositories with AC.txt under 2 UniProt version\n<YYYY_MM_enyo>: UniProt version used by ENYO e.g. 2019_01\n<YYYY_MM_nextprot>: UniProt version used by neXtProt e.g. 2020_01\n<vp_dict_file>: file where the viral protein dictionnary object is stored\t<hp_dict_file>: file where the human protein dictionnary object is stored\n' % sys.argv[0])

if __name__ == "__main__":
    update_ac_uniprot_information(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    #./update_ac_uniprot_information.py vinland_flatfile.tsv results/2019_01/ results/2019_09/
