#AUTHOR: RATTINA Vimel - 2019/11/28 - SIB & Enyo Pharma

import sys #I/O files
import os, errno #create folder
import requests #wget
import contextlib #close the url
import csv #parse tsv
from joblib import Parallel, delayed #parallelize
import multiprocessing #parallelize
import logging #log file

logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', filename='main.log', level=logging.DEBUG)

logging.getLogger("requests").setLevel(logging.WARNING) #make silent requests.get prints
logging.getLogger("urllib3").setLevel(logging.WARNING)

##The aim of this script is to download AC.txt (and AC.tab) files according to a specific version of UniProt.

## INPUTS: ##
##-The flatten PPi file from ENYO (output of flatten_enyo_json.py),
##-An YYYY_MM UniProt version e.g. 2019_01 (version of ENYO or neXtProt),
##-A output folder name,
##-And a "True" or "False" string to know if there is a need to download the AC.tab file or no.

## OUTPUTS: ##
##-Creation of the output folder,
##-A history folder (if "True" specified) with AC.tab files,
##-And a YYYY_MM folder with AC.txt files with the specified UniProt version.

###################
# FUNCTIONS #######
###################

### Function to download every AC.tab in UniProt
def wget_historic(AC, output_dir):
    AC_history = "https://www.uniprot.org/uniprot/"+AC+".tab?version=*"
    
    ##time.sleep(10) #sleep could be useful to let time before to request URL but slow down the process
    tab_content = requests.get(AC_history)
    tab_output = output_dir+"/history/"+AC+".tab"
    
    try:
        tab_content.raise_for_status()
        
        with open(tab_output, 'w') as ftab:
            ftab.write(tab_content.text)
            
        return AC+" downloaded"
        
    except requests.exceptions.RequestException as e:
        ##If failed once retry
        tab_content = requests.get(AC_history)
        with open(tab_output, 'w') as ftab:
            ftab.write(tab_content.text)
            return AC+" downloaded"
        ##If not correct again
        logging.warning(AC_history+"\t404 not found, sorry this page was not found in UniProt")
        return AC+" failed"


### Function to download every AC.txt in UniProt with a specific version
def wget_ac(AC, output_dir, YYYY_MM):
    #print AC
    output_dir_version = output_dir+"/"+YYYY_MM

    AC_tab = output_dir+"/history/"+AC+".tab"
    ##If historic file was downloaded
    if os.path.isfile(AC_tab):
        with open(AC_tab, "r") as tab_file:
            header = tab_file.readline()
            correct_version = False
            entry_v = ""
            while ( correct_version == False ):
                history_line = tab_file.readline().rstrip().split("\t")
                history_date = history_line[4]
                
                if (history_date <= YYYY_MM):
                    correct_version = True 
                    entry_v = history_line[0]

                    ##Sometimes the first line version is 0 because the AC get replaced by a new AC in the latest version
                    ##If YYYY_MM already inferior or equal, do not download the AC_txt because will create an empty file
                    ##e.g. https://www.uniprot.org/uniprot/P30443.tab?version=* for 2019_09
                    if entry_v == "0":
                        if (len(history_line) == 8):
                            replacedby = history_line[7]
                            logging.warning(AC+" is replaced by "+replacedby+" in "+YYYY_MM)
                            return None ##do not try to download a outdate version, get out from the function
                        else: ##not replaced but deleted because obsolete
                            logging.warning(AC+" is obsolete and deleted in "+YYYY_MM)
                            return None
        #print AC+"\t"+str(YYYYhistory)+"\t"+str(MMhistory)+"\t"+str(entry_v)+"\tvs\t"+str(YYYYhistory)+"\t"+str(MMhistory)
   
        ##Download corresponding txt version
        AC_txt="https://www.uniprot.org/uniprot/"+AC+".txt?version="+entry_v

        txt_content = requests.get(AC_txt) ##timeout = None option is a bad solution can take a long time
        try:
            txt_content.raise_for_status()
            with open(output_dir_version+"/"+AC+".txt", 'w') as ftxt:
                ftxt.write(txt_content.text)
        except requests.exceptions.RequestException as e:
            ##if failed once retry just in case
            txt_content = requests.get(AC_txt)
            with open(output_dir_version+"/"+AC+".txt", 'w') as ftxt:
                ftxt.write(txt_content.text)
                return None
            ##if not correct again the page probably does not exist
            logging.warning(AC_txt+"\t404 not found, sorry this page was not found in UniProt")


### Function highlighting the empty files or with only the header to be download again
def print_download_failure(output_dir, file_extension):
    downloaded_files = os.listdir(output_dir)
    AC_failed = []
    for i in downloaded_files:
        n_line = sum(1 for line in open(output_dir+"/"+i)) #without with
        if (n_line == 1) or (n_line == 0): #only header get downloaded or empty
            AC_ID = i.split(file_extension)[0]
            AC_failed.append(AC_ID)
    return AC_failed 
    

def download_ac_in_uniprot_versions(flatten_ppi_tsv, YYYY_MM, output_dir, gethistoric="True"):

    logging.info("download_ac_in_uniprot_versions function starts to run")

    AC_list=[]
    with open(flatten_ppi_tsv) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')

        for line in reader:
            AC_list.append(line['interactor1_accession']) 
            AC_list.append(line['interactor2_accession']) ##reminder: if VH PPi type, interactor2 is the viral protein

    AC_list = list(set(AC_list)) ##remove the duplicated AC
    #print "download order:\t"+str(AC_list)+"\n"
    num_cores = int((multiprocessing.cpu_count() - 3)) #1

    status=[]
    if gethistoric == "True":
        logging.info("AC.tab (history) download process starts")
        history_folder=output_dir+"/history"
        if not os.path.exists(history_folder):
            os.makedirs(history_folder)
        else:
            logging.warning(history_folder+" already exists, files inside will be overwritted")

        status += Parallel(n_jobs=num_cores)(delayed(wget_historic)(i, output_dir) for i in AC_list)
        
        for i in status:
            if " failed" in i:
                print(i+" was not downloaded due to some troubles and will not be integrated in the study")
                ##can put into the log if needed ##here the AC does not exist in UniProt, fake AC
        logging.info("AC.tab (history) download process ends")

        ##check if downloaded correctly
        logging.info("AC.tab (history) re-download process if failed starts")
        AC_redownload = print_download_failure(history_folder, ".tab")
        logging.warning("These AC failed and are in redownloading step:\t"+str(AC_redownload))
        status += Parallel(n_jobs=1)(delayed(wget_historic)(i, output_dir) for i in AC_redownload)
        ##check if redownloaded correctly
        AC_redownload = print_download_failure(history_folder, ".tab")
        logging.warning("These AC failed again:\t"+str(AC_redownload))
        logging.info("AC.tab (history) re-download process if failed ends")

    logging.info("AC.txt download process starts for "+YYYY_MM)
    output_dir_version = output_dir+"/"+YYYY_MM
    if not os.path.exists(output_dir_version):
        os.makedirs(output_dir_version)
    else:
        logging.warning(output_dir_version+" already exists, files inside will be overwritted")

    ACjob = Parallel(n_jobs=num_cores)(delayed(wget_ac)(i, output_dir, YYYY_MM) for i in AC_list)
    logging.info("AC.txt download process ends for "+YYYY_MM)

    ##check if downloaded correctly
    logging.info("AC.txt re-download process if failed starts for "+YYYY_MM)
    AC_redownload = print_download_failure(output_dir_version, ".txt")
    logging.warning("These AC failed and are in redownloading step:\t"+str(AC_redownload))
    ACjob_retry = Parallel(n_jobs=1)(delayed(wget_ac)(i, output_dir, YYYY_MM) for i in AC_redownload)
    ##check if redownloaded correctly
    AC_redownload = print_download_failure(output_dir_version, ".txt")
    logging.warning("These AC failed again:\t"+str(AC_redownload))
    logging.info("AC.txt re-download process if failed ends for "+YYYY_MM)

####################
# MAIN #############
####################

if len(sys.argv) < 5:
    sys.exit('Usage: %s <flat_file> <MM_YYYY> <output_dir> <gethistoric>\n<flat_file>: flatten json file from ENYO Pharma\n<YYYY_MM>: digit(s) for the year_month of the previous UniProt version e.g. 2019_01 for january 2019 UniProt version\n<output_dir>: folder where AC.txt will be stored\n<gethistoric>: for the first use put true to download the historic otherwise it will not work, for the second use (accessory the second date) historic is already downloaded so use FALSE to avoid waste of time' % sys.argv[0])

if __name__ == "__main__":        
    download_ac_in_uniprot_versions(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
#./download_ac_in_uniprot_versions.py EXAMPLE.tsv 2019_01 results True
#./download_ac_in_uniprot_versions.py EXAMPLE.tsv 2019_09 results False
