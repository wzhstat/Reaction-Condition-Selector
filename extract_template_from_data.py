import pandas as pd
import re
import numpy as np
from tqdm import tqdm
import template_extractor
from rxnmapper import BatchedMapper
from func_timeout import func_timeout, FunctionTimedOut
from joblib import Parallel, delayed
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdChemReactions import RemoveMappingNumbersFromReactions
from rdkit import Chem
import argparse

def remove_chirality(reaction):
    reagent, product = reaction.split('>>')
    reagent = reagent.split('.')
    product = product.split('.')
    reagent = [Chem.MolFromSmiles(i) for i in reagent]
    product = [Chem.MolFromSmiles(i) for i in product]
    reagent = [Chem.MolToSmiles(i, isomericSmiles=False) for i in reagent]
    product = [Chem.MolToSmiles(i, isomericSmiles=False) for i in product]
    reagent = '.'.join(reagent)
    product = '.'.join(product)
    smiles = reagent + '>>' + product
    return smiles

def remove_mapping(reaction):
    try:
        reaction = rdChemReactions.ReactionFromSmarts(reaction)
        RemoveMappingNumbersFromReactions(reaction)    
        reaction = rdChemReactions.ReactionToSmiles(reaction)
        return reaction
    except:
        return ''

def _extract(reaction,radius=1):
    try:
        return template_extractor.extract_from_reaction(reaction,radius)
    except KeyboardInterrupt:
        print('Interrupted')
        raise KeyboardInterrupt
    except Exception:
        return {'reaction_id': reaction['_id']}

def extract(reaction,radius):
    try:
        return func_timeout(20, _extract, args=(reaction,radius))
    except FunctionTimedOut:
        print('Timeout')
        return {'reaction_id': reaction['_id']}

def can_parse(reaction):
    try:
        reagenr, product = reaction.split('>>')
        if Chem.MolFromSmiles(reagenr) and Chem.MolFromSmiles(product):
            return True
        else:
            return False
    except:
        return False

def get_data(args):
    data = pd.read_csv(args.data_path)
    can_parse_list = Parallel(n_jobs=-1)(delayed(can_parse)(reaction) for reaction in data[args.reaction_smiles_column])
    rxn_mapper = BatchedMapper(batch_size=32)
    # skip the reactions that cannot be parsed and don't map the reactions(but don't remove them,fill the mapped reactions with the original reactions)
    data = data[can_parse_list] 
    data[args.reaction_smiles_column] = [remove_chirality(reaction) for reaction in data[args.reaction_smiles_column]]
    print('Mapping reactions...')
    data['Mapped_Reaction'] = list(rxn_mapper.map_reactions(list(data[args.reaction_smiles_column])))
    print('Done')
    split_smiles = data['Mapped_Reaction'].str.split('>>',expand=True)
    data['reactants'] = split_smiles[0]
    data['products'] = split_smiles[1]
    data['_id'] = data[args.id_column]
    reactions = data[['_id', 'reactants', 'products']].to_dict('records')
    print('Extracting templates...')
    templates_r0 = Parallel(n_jobs=-1, verbose=4)(delayed(extract)(reaction,0) for reaction in reactions)
    templates_r1 = Parallel(n_jobs=-1, verbose=4)(delayed(extract)(reaction,1) for reaction in reactions)
    print('Done')


    # Add the templates to the DataFrame
    templates_list0 = []
    for template in templates_r0:
        try:
            templates_list0.append(template['reaction_smarts'])
        except:
            templates_list0.append('')
        # Add the templates to the DataFrame
    templates_list1 = []
    for template in templates_r1:
        try:
            templates_list1.append(template['reaction_smarts'])
        except:
            templates_list1.append('')
    data['template_r0'] = templates_list0
    data['template_r1'] = templates_list1
    templates_r_1 = [remove_mapping(tpl) for tpl in data['template_r0']]
    data['template_r0*'] = templates_r_1

    # drop if the template is empty
    data = data[data['template_r0'] != '']
    data = data[data['template_r1'] != '']
    data = data[data['template_r0*'] != '']

    data.to_csv(args.out_path, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', type=str, default='test_template_extractor.csv')
    parser.add_argument('--reaction_smiles_column', type=str, default='reaction')
    parser.add_argument('--id_column', type=str, default='_id')
    parser.add_argument('--out_path', type=str, default='out.csv')
    args = parser.parse_args()
    get_data(args)