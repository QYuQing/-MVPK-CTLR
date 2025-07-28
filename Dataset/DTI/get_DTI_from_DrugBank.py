import os
import xml.etree.ElementTree as ET
from tqdm import tqdm

############################################# THE Best ONE #############################################################
def get_DTI_coordinates(drug_entry, coordinate_list):
    drugbank_ID = drug_entry.find('drugbank-id').text
    ## TARGETS
    if drug_entry.findall('.//targets')[0]: 
        get_target(drug_entry, coordinate_list, drugbank_ID)
    # ENZYMES
    if drug_entry.findall('.//enzymes')[0]: 
        get_enzyme(drug_entry, coordinate_list, drugbank_ID)
    ## TRANSPORTERS
    if drug_entry.findall('.//transporters')[0]:
        get_transporter(drug_entry, coordinate_list, drugbank_ID)
    ## CARRIERS
    if drug_entry.findall('.//carriers')[0]:
        get_carrier(drug_entry, coordinate_list, drugbank_ID)
    return coordinate_list
##
def get_target(drug_entry, coordinate_list, drugbank_ID):
    target = None
    for tgt in  drug_entry.findall('.//targets')[0]: # iterate in targets
        for i in tgt.findall('.//polypeptide'): # for searching external-identifiers
            ext = i.findall('.//external-identifiers')[0]
            for id in list(ext):
                if list(id)[0].text == "UniProtKB": # can be done with get
                    target = list(id)[1].text
                    if target:
                        coordinate = (drugbank_ID, target)
                        coordinate_list.append(coordinate)
## 
def get_enzyme(drug_entry, coordinate_list, drugbank_ID):
    list_enzymes = drug_entry.findall('.//enzymes')[0]
    enzyme = None
    for enz in list(list_enzymes):
        if enz.findall('.//uniprot-id'):
            enzyme = enz.findall('.//uniprot-id')[0].text
        elif enz.findall('.//polypeptide'):
            pol = enz.findall('.//polypeptide')[0]
            if pol.findall('.//external-identifiers')[0]:
                exts = pol.findall('.//external-identifiers')[0]
                for id in list(exts):
                    if list(id)[0].text == "UniProtKB": # can be done with get
                        enzyme = list(id)[1].text
        if enzyme:
            coordinate = (drugbank_ID, enzyme)
            coordinate_list.append(coordinate)
##
def get_transporter(drug_entry, coordinate_list, drugbank_ID):
    list_transporters = drug_entry.findall('.//transporters')[0]
    transporter = None
    for trans in list(list_transporters):
        if trans.find(('.//polypeptide')):
            transporter = trans.find(('.//polypeptide')).get('id')
        if transporter:
            coordinate = (drugbank_ID, transporter)
            coordinate_list.append(coordinate)
##
def get_carrier(drug_entry, coordinate_list, drugbank_ID):
    list_carriers = drug_entry.findall('.//carriers')[0]
    carrier = None
    for car in list(list_carriers):
        if car.find(('.//polypeptide')):
            carrier = car.find(('.//polypeptide')).get('id')
        if carrier:
            coordinate = (drugbank_ID, carrier)
            coordinate_list.append(coordinate)




# read xml
tree = ET.parse(r'F:\DataBase\DrugBank\full database.xml')
# tree = ET.parse(r'F:\DataBase\DrugBank\DrugBankParser-master\main\test50.xml')
root = tree.getroot()

# execute the function
coordinate_list = []
for drug_entry in tqdm(root):
    coordinate_list = get_DTI_coordinates(drug_entry, coordinate_list)

print(f'Coordinate list len: {len(coordinate_list)}, from those {len(set(coordinate_list))} are unique')

db_path = os.path.join(r'F:\DataBase\DrugBank\DTI\DrugBank_DTIs.tsv')
with open(db_path, 'w') as f:
    _ = f.write('# DrugBank ID\tUniprot ID\n')
    for item in coordinate_list:
        _ = f.write("%s\t%s\n" % item)
