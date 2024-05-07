import csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
import numpy as np
import pandas as pd
import sys
from joblib import Parallel, delayed
import csv
import os
import argparse
import gzip
csv.field_size_limit(500 * 1024 * 1024)

def get_keys(path:str):
    '''
    open csv file and get all cat,solv,reag keys
    args:
        path: csv file path
    '''
    
    with open('%s/all_cat_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        cat_list = [row['cat'] for row in reader]
    
    with open('%s/all_solv_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        solv_list = [row['solv'] for row in reader]

    with open('%s/all_reag_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        reag_list = [row['reag'] for row in reader]
    return cat_list,solv_list,reag_list 



def remove_reagent(Smarts:str):
    '''
    remove reagent from reaction
    args:
        Smarts: reaction smarts
    '''
    rxn = Smarts.split('>')
    reactant = rxn[0]
    product = rxn[2]
    reactant = reactant.split('.')
    product = product.split('.')
    outreactant = filter(lambda i: 'C:' in i or 'CH:' in i or 'CH2:' in i or 'CH3:' in i or 'N:' in i or 'NH:' in i or 'NH2:' in i or 'NH3:' in i or 'F:' in i or 'Cl:' in i or 'Br:' in i or 'I:' in i, reactant)
    rout = '.'.join(outreactant)
    outproduct = filter(lambda i: 'C:' in i or 'CH:' in i or 'CH2:' in i or 'CH3:' in i or 'N:' in i or 'NH:' in i or 'NH2:' in i or 'NH3:' in i or 'F:' in i or 'Cl:' in i or 'Br:' in i or 'I:' in i, product)
    pout = '.'.join(outproduct)
    out = rout + '>>' + pout
    return out


def get_index(condition,list):
    try:
        if str(condition) == 'nan':
            condition = "None"
        return list.index(condition)
    except Exception as e:
        return 0

def Extraction_MPNN_data(args,target_list: list, data: pd.DataFrame):
    '''
    open csv file and get all data for MPNN
    args:
        path: csv file path
        file_name: csv file name
        target: target name
        target_list: target list
        data: csv file data
        condition: condition that is used for MPNN
        cat_list: cat list
        solv_list: solv list
        reag_list: reag list
    '''
    #in_lists = Parallel(n_jobs=-1, verbose=4)(delayed(in_list)(condition,target_list) for condition in list(data[args.target]))
    #data = data[in_lists]
    MLP_all_data = pd.DataFrame()
    rxnsmile = Parallel(n_jobs=-1, verbose=4)(delayed(remove_reagent)(reaction) for reaction in list(data['reaction']))
    target_index = Parallel(n_jobs=-1, verbose=4)(delayed(get_index)(condition,target_list) for condition in list(data[args.target]))
    MLP_all_data['reaction'] = rxnsmile
    MLP_all_data['target'] = target_index
    return MLP_all_data


def save_csv(args,out_data):
    '''
    save data to csv file
    args:
        path: save path
        condition: condition that is used for MPNN
        MLP_all_data: all data 
        target: target name
        N: whether to use N
    '''
    data_name = "GCN_%s.csv"%args.data_name
    if os.path.exists(args.save_path):
        pass
    else:
        os.mkdir(args.save_path)
    path = '%s/GCN_data'%(args.save_path)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)
    out_data.to_csv('%s/%s'%(path,data_name),index=False)
    print('save to %s/%s'%(path,data_name))

def get_MPNN_data(args):
    '''
    get data for GCNN
    args:
        path: csv file path
        file_name: csv file name
        target: target name
    '''
    data = pd.read_csv('%s/%s.csv'%(args.data_path,args.data_name))
    cat_list,solv_list,reag_list = get_keys(args.key_path)
    if args.target in ['cat']:
        target_list = cat_list
    elif args.target in ['solv0','solv1']:
        target_list =solv_list
    else:
        target_list = reag_list
    MLP_all_data= Extraction_MPNN_data(args,target_list,data)
    return MLP_all_data

if __name__ == '__main__':
    pass

    