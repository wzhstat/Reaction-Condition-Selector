from rdchiral import template_extractor
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdChemReactions import RemoveMappingNumbersFromReactions
import hashlib
import pandas as pd

def get_tem(reaction,n = 0):
    data = {}
    data['_id'] = n
    split_smiles =reaction.split('>')
    data['reactants'] = split_smiles[0]
    data['spectators'] = split_smiles[1]
    data['products'] = split_smiles[2] 
    
    try:
        template = template_extractor.extract_from_reaction(data)
        reaction_smarts = template['reaction_smarts']
        if reaction_smarts is not None:
            reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
            RemoveMappingNumbersFromReactions(reaction)
            reaction = rdChemReactions.ReactionToSmiles(reaction)
        return reaction
    except Exception:
        return 'else'


def binary_search(key,hash_key):
    low = 0
    high = len(hash_key)-1
    while low <= high:
        mid = (low+high)//2
        guess = hash_key[mid]
        if guess == key:
            return mid
        else:
            if guess > key:
                high = mid -1
            if guess < key:
                low = mid + 1
    return 0

def popularity_module0(input_data = None):
    return [0,1,2,3,4,5,6,7,8,9],[0,1,2,3,4,5,6,7,8,9]



def popularity_module1(input_data,popularity_data = pd.read_csv('data/popular_data.csv')):
    reaction = input_data
    hash_key = list(popularity_data['hash_id'])
    tem = get_tem(reaction)
    hash_id = hashlib.md5(tem.encode()).hexdigest()
    n = binary_search(hash_id,hash_key)
    cat_out = {}
    solv_out = {}
    for i in range(10):
        cat_out['cat_top-%d'%(i+1)] = popularity_data['cat_top-%d'%(i+1)][n]
        solv_out['solv_top-%d'%(i+1)] = popularity_data['solv_top-%d'%(i+1)][n]
    return cat_out,solv_out


