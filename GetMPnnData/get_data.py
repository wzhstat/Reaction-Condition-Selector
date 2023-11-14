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
import gzip
csv.field_size_limit(500 * 1024 * 1024)

def open_csv(path:str):
    '''
    open csv file and get all cat,solv,reag
    args:
        path: csv file path
    '''
    with open('%s/all_cat_withoutN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            cat_list = list(classes.keys())
    
    with open('%s/all_cat_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            cat_list_N = list(classes.keys())

    with open('%s/all_solv_withoutN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            solv_list = list(classes.keys())
    
    with open('%s/all_solv_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            solv_list_N = list(classes.keys())

    with open('%s/all_reag_withoutN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            reag_list = list(classes.keys())

    with open('%s/all_reag_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            reag_list_N = list(classes.keys())
    return cat_list,cat_list_N,solv_list,solv_list_N,reag_list,reag_list_N 



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
    outreactant = filter(lambda x: 'C:' in x or 'CH:' in x or 'CH2:' in x or 'CH3:' in x, reactant)
    rout = '.'.join(outreactant)
    outproduct = filter(lambda x: 'C:' in x or 'CH:' in x or 'CH2:' in x or 'CH3:' in x, product)
    pout = '.'.join(outproduct)
    out = rout + '>>' + pout
    return out

def in_list(condition,list):
    if str(condition) == 'nan':
        condition = "None"
    if condition in list:
        return True
    else:
        return False
def get_index(condition,list):
    if str(condition) == 'nan':
        condition = "None"
    return list.index(condition)

def get_condition_fp(condition_siles:str):
    '''
    get condition fp
    args:
        condition_siles: condition smiles
    '''
    try:
        m = Chem.MolFromSmiles(condition_siles)
        fp = np.array(AllChem.GetMorganFingerprintAsBitVect(m, useChirality=True, radius=2, nBits = 512))
        return fp
    except Exception as e:
        return np.zeros(512)

def get_condition_one_hot(condition,list):
    '''
    get condition one-hot
    args:
        condition: condition
        list: condition list
    '''
    get_one_hot = np.zeros(len(list))
    if str(condition) == 'nan':
        condition = "None"
    if condition in list:
        one_hot_index = list.index(condition)
        get_one_hot[one_hot_index] = 1
    return get_one_hot

def open_MPNN_data(path: str,file_name: str,target: str,target_list: list, data: pd.DataFrame,conditions: list, module: str,):
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
    in_lists = Parallel(n_jobs=-1, verbose=4)(delayed(in_list)(condition,target_list) for condition in list(data[target]))
    data = data[in_lists]
    MLP_all_data = pd.DataFrame()
    rxnsmile = Parallel(n_jobs=-1, verbose=4)(delayed(remove_reagent)(reaction) for reaction in list(data['reaction']))
    target_index = Parallel(n_jobs=-1, verbose=4)(delayed(get_index)(condition,target_list) for condition in list(data[target]))
    MLP_all_data['smarts'] = rxnsmile
    MLP_all_data['target'] = target_index
    all_condition_list = []
    if conditions != None:
        _,cat_list,_,solv_list,_,reag_list = open_csv(path)
        if module == 'one-hot':
            all_condition_list = []
            for condition in conditions:
                if condition in ['cat','solv']:
                    condition_lists = eval('%s_list'%condition)
                else:
                    condition_lists = reag_list
                condition_list = Parallel(n_jobs=-1, verbose=4)(delayed(get_condition_one_hot)(condition,condition_lists) for condition in list(data[condition]))
                all_condition_list.append(condition_list)
            all_condition_list = [[clist[i] for clist in all_condition_list] for i in range(len(all_condition_list[0]))]
            all_condition_list = [np.concatenate(clist) for clist in all_condition_list]
        elif module == 'fp':
            all_condition_list = []
            for condition in conditions:
                condition_list = Parallel(n_jobs=-1, verbose=4)(delayed(get_condition_fp)(condition) for condition in list(data[condition]))
                all_condition_list.append(condition_list)
            all_condition_list = [[clist[i] for clist in all_condition_list] for i in range(len(all_condition_list[0]))]
            all_condition_list = [np.concatenate(clist) for clist in all_condition_list]
    return MLP_all_data,all_condition_list



def save_MPNN_csv(path,conditions,features,MLP_all_data,module,target,N):
    '''
    save data to csv file
    args:
        path: save path
        condition: condition that is used for MPNN
        MLP_all_data: all data 
        target: target name
        N: whether to use N
    '''
    header = ['smarts','target']
    data_name = "GCN_data.csv"
    if N:
        path = '%s/GCN_%s_data_withN'%(path,target)
    else:
        path = '%s/GCN_%s_data_withoutN'%(path,target)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)
    MLP_all_data.to_csv(path+"/%s"%(data_name), mode='a', index=False,header=header)
    if conditions != None:
        features = np.array(features)
        np.savez_compressed(path+"/%s.npz"%('+'.join(conditions)+module), features=features)

def get_GCNN_data(path: str, file_name: str, target: str, conditons: list, condition_moduel: str = '', N: bool = False):
    '''
    get data for GCNN
    args:
        path: csv file path
        file_name: csv file name
        target: target name
        conditons: condition that is used for MPNN
        N: whether to use N
    '''
    data = pd.read_csv('%s/%s.csv'%(path,file_name))
    cat_list,cat_list_N,solv_list,solv_list_N,reag_list,reag_list_N = open_csv(path)
    if N:
        if target in ['cat','solv']:
            target_list = eval('%s_list_N'%target)
        else:
            target_list = reag_list_N
    else:
        if target in ['cat','solv']:
            target_list = eval('%s_list'%target)
        else:
            target_list = reag_list
    MLP_all_data,features= open_MPNN_data(path,file_name,target,target_list,data,conditons,condition_moduel)
    save_MPNN_csv(path,conditons,features,MLP_all_data,condition_moduel,target,N)

if __name__ == '__main__':
    get_GCNN_data('./data','1976-2016_5+','cat',None,N=True)
    get_GCNN_data('./data','1976-2016_5+','cat',None)
    get_GCNN_data('./data','1976-2016_5+','solv',['cat'],'fp',N=True)
    get_GCNN_data('./data','1976-2016_5+','solv',['cat'],'fp')
    get_GCNN_data('./data','1976-2016_5+','reag0',['cat','solv'],'fp',N=True)
    get_GCNN_data('./data','1976-2016_5+','reag0',['cat','solv'],'fp')
    get_GCNN_data('./data','1976-2016_5+','solv',['cat'],'one-hot',N=True)
    get_GCNN_data('./data','1976-2016_5+','solv',['cat'],'one-hot')
    get_GCNN_data('./data','1976-2016_5+','reag0',['cat'],'one-hot',N=True)
    get_GCNN_data('./data','1976-2016_5+','reag0',['cat'],'one-hot')

    
