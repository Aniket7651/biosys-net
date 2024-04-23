"""
Created on Aug 20 11:21:32 2022
This is the file provided by pyCrossbill old pack, for reteriving the data from external biological database site like, 
DrugBank by programetically 
@author: ANIKET YADAV
"""


import os
import re
from bs4 import BeautifulSoup
import requests
from urllib.request import urlopen
from urllib.parse import quote


exception = 'ERR: reterive failed:\n1. check your internet connection first\n2. wrong access database or mistakes!'

# protein_url = 'https://www.ncbi.nlm.nih.gov/protein/?term=cry'
# nucleotide_url = 'https://www.ncbi.nlm.nih.gov/nuccore/?term=cry'
# gene_url = https://www.ncbi.nlm.nih.gov/gene/?term=cry
# out = urlopen(url).read().decode('utf-8')

def get_ID(query, database_name='protein'):
    '''getting unique ID of the query from NCBI (protein, gene, neucleotide) database
    it's take query(what do you want to search about) and database name, where default database 'protein' is selected'''

    if database_name == 'neucleotide':
        url = f'https://www.ncbi.nlm.nih.gov/nuccore/?term={query}'
        
    else:
        url = f'https://www.ncbi.nlm.nih.gov/{database_name}/?term={query}'

    try:
        soup = BeautifulSoup(requests.get(url).text, 'lxml')
        id_list = []
        for dd in soup.select('.rprtid dd'):
            id_list.append(dd.string)
        return id_list

    except requests.exceptions.ConnectionError:
        return exception


def drug_ID(drug_name):
    '''Search for drug by drug name from drug bank database\nreturns very first drug ID from drug list'''

    try:
        url = f'https://go.drugbank.com/unearth/q?query={drug_name}&button=&searcher=drugs'
        soup = BeautifulSoup(requests.get(url).text, 'lxml')
        id_list = []
        for dd in soup.select('.hit-link a') or soup.find('dd', {'class': 'col-xl-4 col-md-9 col-sm-8'}).get_text():
            id_list.append(re.findall("DB\d\d\d\d\d", str(dd)))
        return id_list[0][0]
    except TypeError:
        return 'probabilly wrong or spell mistake'
    except requests.exceptions.ConnectionError:
        return exception


def DID(drug_name):
    '''Search for drug by drug name from drug bank database\nreturns very first drug ID from drug list
    \nupdated method of 'drug_ID' function in some times drug_ID method doesn't work, in this case you can use this method'''

    try:
        url = f'https://go.drugbank.com/unearth/q?query={drug_name}&button=&searcher=drugs'
        soup = BeautifulSoup(requests.get(url).text, 'lxml')
        id_list = []
        for dd in soup.find('dd', {'class': 'col-xl-4 col-md-9 col-sm-8'}):
            id_list.append(dd.get_text())
        return id_list[0]
    except TypeError:
        return 'probabilly wrong or spell mistake'
    except requests.exceptions.ConnectionError:
        return exception


class drug_info():

    def __init__(self, drugBank_ID):
        self.id = drugBank_ID

    def drug_detail(self):
        try:
            url = 'https://go.drugbank.com/drugs/%s' %self.id
            soup = BeautifulSoup(requests.get(url).text, 'lxml')
            (key_lis, value_lis) = ([], [])
            for dt, dd in zip(soup.select('.card-content dl dt'), soup.select('.card-content dl dd')):
                key_lis.append(dt.string)
                value_lis.append(dd.string)
    
            data = {key_lis[i]: value_lis[i] for i in range(len(key_lis))}
            return data
        except requests.exceptions.ConnectionError:
            return exception
    
    def smile(self):
        data = self.drug_detail()
        return data['SMILES']

    def CAS_number(self):
        data = self.drug_detail()
        return data['CAS number']

    def IUPAC_name(self):
        data = self.drug_detail()
        return data['IUPAC Name']

    def InchI_key(self):
        data = self.drug_detail()
        return data['InChI Key']

    def PubChem_id(self):
        data = self.drug_detail()
        return data['PubChem Compound']


class targets_info(drug_info):

    def __init__(self, drug_id):
        self.id = drug_id
        # self.dinfo = drug_info(self).drug_detail()

    def drug_targets(self):
        try: 
            url = f'https://go.drugbank.com/drugs/{self.id}'
            soup = BeautifulSoup(requests.get(url).text, 'lxml')
            target_str = soup.find('table',{'class':'table table-sm responsive-table'}).get_text(',')
            list_target = target_str.split(',')
            target = []
            for i in range(3, len(list_target), 4):
                target.append(list_target[i:i+4])
            return target
        except requests.exceptions.ConnectionError:
            return exception

    def get_information(self):
        try:
            url = f'https://go.drugbank.com/drugs/{self.id}'
            soup = BeautifulSoup(requests.get(url).text, 'lxml')
            (dt_lis, dd_lis) = ([], [])
            for dt, dd in zip(soup.select('.bond dl dt'), soup.select('.bond dl dd')):     # strong a
                (dt_lis.append(dt.string), dd_lis.append(dd.string))
            
            data = []
            for i, j in zip(dt_lis, dd_lis):
                data.append(f'{i}: {j}')
            return data
        except requests.exceptions.ConnectionError:
            return exception

    def GetUniprotID(self):
        get_info = self.get_information()
        for id in get_info:
            if('Uniprot ID' in id):
                return id.split(':')[1][1:] # return last or second value from second index by [1][1:]

    def target_file(self, tsv=False):
        dinfo = drug_info(self.id).drug_detail()
        info = self.get_information()
        if(tsv == True):
            with open(f'{self.id}_target_info.tsv', 'a') as target_file:
                target_file.writelines('Type\t\tOrganism\t\tPharmacological action\t\tActions\t\tGeneral Function\t\tSpecific Function\t\tGene Name\t\tUniprot ID\t\tUniprot Name\t\tMolecular Weight\n')
                count = 0
                for l in range(len(info)):  
                    if(count == 9):
                        count = 0 
                        target_file.write('\n')
                    else:
                        target_file.write("%s\t\t\t"%info[l].split(':')[1][1:])
                        count += 1
            return os.path.abspath(f'{self.id}_target_info.tsv')

        else:
            with open(f'{self.id}_target_info.trg', 'a') as target_file:
                co = 0
                target_file.write('$INFO>>{}//{}//{}// {} || {} || {};\n'.format(self.id, dinfo['Generic Name'], dinfo['IUPAC Name'], dinfo['Super Class'], dinfo['Class'], dinfo['Sub Class']))
                for i in info[1:]:
                    if co < 9:
                        target_file.write(i)
                        target_file.write('\n')
                        co += 1
                    else:
                        target_file.write('\n\n')
                        co = 0
                target_file.write('<<TERMINATED$')
            return os.path.abspath(f'{self.id}_target_info.trg')


def getSimilarTargets(cid):
    try:
        url = 'https://go.drugbank.com/polypeptides/%s' %cid
        soup = BeautifulSoup(requests.get(url).text, 'lxml')
        table = soup.findAll('table')  # find all data with HTML tag 'table'
        html = table[table.__len__()-1] # accept last most element in the list
        # header and fdata(full data) stores content which exist in tag 'tr'(table row) and 'td'(table data)
        header = [th.text for th in html.findAll('th')] 
        fdata = [td.text for td in html.findAll('td')]
        # exactly, fdata store all the row of HTML table in a single list it is nessesery to store it by rows 
        # doing this by making a loop and append numbers of data up to its header length in the sublist and store sublist in 'data' list
        data = []
        for i in range(0, len(fdata), len(header)):
            data.append((fdata[i:i+len(header)]))   # select up to len of header
        return header, data # return data as 2D list and header as list
    except requests.exceptions.ConnectionError:
        return exception
    # soup = BeautifulSoup(requests.get(url).text, 'lxml')
    # seq = soup.find('pre', {'class':'col-xl-10 col-md-9 col-sm-8'})
    # return seq


# Convert some drug data to table form return two output header of table and 2D list of info data
# there are the HTML tables in a single scraping so, tableLoader loads individual HTML table
def DrugTrialsTable(drugID, tableLoader):
    try:
        soup = BeautifulSoup(requests.get('https://go.drugbank.com/drugs/%s'%drugID).text, 'lxml')
        htm = soup.findAll('table')[tableLoader]  # find all data with HTML tag 'table'
        # header and fdata(full data) stores content which exist in tag 'tr'(table row) and 'td'(table data)
        header = [th.text for th in htm.findAll('th')] 
        fdata = [td.text for td in htm.findAll('td')]
        # exactly, fdata store all the row of HTML table in a single list it is nessesery to store it by rows 
        # doing this by making a loop and append numbers of data up to its header length in the sublist and store sublist in 'data' list
        data = []
        for i in range(0, len(fdata), len(header)):
            data.append((fdata[i:i+len(header)]))   # select up to len of header
        return header, data # return data as 2D list and header as list
    except requests.exceptions.ConnectionError:
        return exception


def toCSV(array, header, outfile):
    with open(outfile, 'w') as f:
        for h in header:
            f.write('%s,'%h)
        for i in range(len(array)):
            for j in range(len(array[0])):
                f.write('%s,'%array[i][j])
            f.write('\n')


############################# END OF THE PROGRAM AND COMMANDS CALL ######################################

import sys
task = sys.argv[1] # drugb-[*task]
did = sys.argv[2] # drug id

# ~$ drugb-target DB00291 -tsv OR -csv
if(task == "drugb-similardrug"):
    smilar = getSimilarTargets(did)
    toCSV(smilar[1], smilar[0], sys.argv[3]) # output file name

elif(task == "drugb-target" and sys.argv[3] == '-tsv'):
    targets_info(did).target_file(tsv=True)
else: targets_info(did).target_file()

# print(get_ID('rac A', database_name='nucleotide')[0])
# drug =  'paractamol' #'morphine'
# data = drug_info('DB00316').drug_detail()
# drug_IDs = DID('paracetamol')
# print(drug_IDs)
# uniprot = targets_info(drug_IDs).GetUniprotID()
# print(getSimilarTargets(uniprot))
# print(DrugTables(drug_IDs, 9))
# print(data['CAS number'])
# smilar = getSimilarTargets('DB00316')
# toCSV(smilar[1], smilar[0], 'outf.csv')
# print(drug_info('DB00316').drug_detail())