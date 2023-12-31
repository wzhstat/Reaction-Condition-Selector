import json
import gzip
import csv
import hashlib
import pandas as pd
from rdkit import Chem
from joblib import Parallel, delayed
from time import time
from rdchiral import template_extractor
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdChemReactions import RemoveMappingNumbersFromReactions


def can_parse(rsmi):
    """
    Check if a SMILES string can be parsed.

    Args:
        rsmi (str): The SMILES string to check.

    Returns:
        bool: True if the SMILES string can be parsed, False otherwise.
    """
    #rea = rdChemReactions.ReactionFromSmarts(str(rsmi))
    #RemoveMappingNumbersFromReactions(rea)
    #rea = rdChemReactions.ReactionToSmiles(rea)
    
    # Get reactants and products
    react, spec, prod = rsmi.split('>')

    # Check if reactants and products can be parsed
    if Chem.MolFromSmiles(react) and Chem.MolFromSmiles(prod):
        return True
    else:
        return False

def extract(reaction):
            try:
                return template_extractor.extract_from_reaction(reaction)
            except KeyboardInterrupt:
                print('Interrupted')
                raise KeyboardInterrupt
            except Exception:
                return {'reaction_id': reaction['_id']}

    
def get_temp(data_name,save_path = "./data"):
    data_path = "%s/%s.csv"%(save_path,data_name)
    save_path1 = "%s/templates1.json.gz"%save_path
    save_path2 = "%s/templates2.json.gz"%save_path
    '''
    This function is used to get the templates from the data.

    Args:
        data_name (str): The name of the data.
        save_path (str): The path to save the template data.

    '''

    # Read data from csv file
    datas = pd.read_csv(data_path, chunksize=10000)

    # Extract templates from reactions
    for data in datas:
        data['ReactionSmiles'] = data['reaction']
        split_smiles =data['ReactionSmiles'].str.split('>', expand=True)
        data['reactants'] = split_smiles[0]
        data['products'] = split_smiles[2]   
        parsable = Parallel(n_jobs=-1, verbose=4)(delayed(can_parse)(rsmi) for rsmi in data['ReactionSmiles'].values)
        data = data[parsable]
        data['_id'] = data['_id']
        reactions = data[['_id', 'reactants', 'products']].to_dict('records')
        templates = Parallel(n_jobs=-1, verbose=4)(delayed(extract)(reaction) for reaction in reactions)
        print('done extract')
        templates = pd.DataFrame(templates)
        templates.to_json("%s/templates.json.gz"%(save_path), orient='records', compression='gzip')


def classif(template):
    _id = template['reaction_id']
    reaction_smarts = template['reaction_smarts']
    if reaction_smarts is not None:
        if reaction_smarts not in adic:
            adic[reaction_smarts] = [_id]
        else:
            adic[reaction_smarts].append(_id)
    
        
def classif_by_temp(data_name,save_path = "./data"):
    get_temp(data_name,save_path)
    temp_path = "%s/templates.json.gz"%save_path
    out_path = "%s/classif_by_temp.csv"%save_path
    '''
    This function is used to classify the data by templates.
    '''
    global adic
    adic = {}
    with gzip.open(temp_path) as f:
        templates = json.load(f)
    for template in templates:
        classif(template)
    print('done classify')

    #print(len(adic))

    bdic = {'else':[]}
    for i in adic:
        if len(adic[i]) >=5:
            bdic[i] = adic[i]
        else:
            for j in adic[i]:
                bdic['else'].append(j)
    keys = bdic.keys()
    with open(out_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile,fieldnames=keys)
        writer.writeheader()
        writer.writerow(bdic)

