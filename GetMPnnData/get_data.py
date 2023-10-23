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
            print(len(cat_list))
    
    with open('%s/all_cat_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            cat_list_N = list(classes.keys())
            print(len(cat_list_N))

    with open('%s/all_solv_withoutN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            solv_list = list(classes.keys())
            print(len(solv_list))
    
    with open('%s/all_solv_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            solv_list_N = list(classes.keys())
            print(len(solv_list_N))

    with open('%s/all_reag_withoutN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            reag_list = list(classes.keys())
            print(len(reag_list))

    with open('%s/all_reag_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            reag_list_N = list(classes.keys())
            print(len(reag_list_N))
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
            


def open_MPNN_data(path: str,file_name: str,target: str,target_list: list, data: pd.DataFrame,condition: list, cat_list,solv_list,reag_list):
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
    data = data
    rxnsmile = Parallel(n_jobs=-1, verbose=4)(delayed(remove_reagent)(reaction) for reaction in list(data['reaction']))
    MLP_all_data = []
    all_condition_list = []
    l = len(rxnsmile)
    for i in range(l):
        if i % (l//100) == 0:
            print("\r", end="")
            print("Progress: {}%: ".format(i/(l//100)), "▋" * int(i/(l//50)), end="")
            sys.stdout.flush()
        dic = {}
        if data[target][i] is np.nan:
            titem = 'None'
        else:
            titem = data[target][i]
        if titem not in target_list:
            continue
        dic['smarts'] = rxnsmile[i]
        dic[target] = target_list.index(titem)
        condition_list = []
        if condition != None:
            for conditions in condition:
                get_one_hot = np.zeros(len(eval('%s_list'%conditions)))
                if data[conditions][i] is np.nan:
                    citem = 'None'
                else:
                    citem = data[conditions][i]
                if citem in eval('%s_list'%conditions):
                    one_hot_index = eval('%s_list'%conditions).index(citem)
                    get_one_hot[one_hot_index] = 1
                condition_list.append(get_one_hot)
            if len(condition_list) == 0:
                continue
            if len(condition_list) > 1:
                ohc = np.concatenate(condition_list)
            else:
                ohc = condition_list[0]
            all_condition_list.append(ohc)
        MLP_all_data.append(dic)
    condition_np = np.array(all_condition_list)
    return MLP_all_data,condition_np

def save_MPNN_csv(path,condition,MLP_all_data,condition_np,target,N):
    '''
    save data to csv file
    args:
        path: save path
        condition: condition that is used for MPNN
        MLP_all_data: all data
        condition_np: condition data
        target: target name
        N: whether to use N
    '''
    if condition == None:
        header = ['smarts',target]
    else:
        header = ['smarts',target]
        condition_header = [','.join(condition)]
    if N:
        path = '%s/GCN_%s_data_withN'%(path,target)
    else:
        path = '%s/GCN_%s_data_withoutN'%(path,target)
    os.mkdir(path)
    with open(path+"/GCN_data.csv",'w',newline='') as f:
        f_csv = csv.DictWriter(f,header,extrasaction='ignore')
        f_csv.writeheader()
        f_csv.writerows(MLP_all_data)
    if condition != None:
        np.save(path+"/condition.npy",condition_np)

def get_GCNN_data(path: str, file_name: str, target: str, conditons: list, N: bool = False):
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
        data,cat_list,solv_list,reag_list = data,cat_list_N,solv_list_N,reag_list_N
    else:
        if target in ['cat','solv']:
            target_list = eval('%s_list'%target)
        else:
            target_list = reag_list
        data,cat_list,solv_list,reag_list = data,cat_list,solv_list,reag_list
    MLP_all_data,condition_np = open_MPNN_data(path,file_name,target,target_list,data,conditons,cat_list,solv_list,reag_list)
    save_MPNN_csv(path,conditons,MLP_all_data,condition_np,target,N)

if __name__ == '__main__':
    get_GCNN_data('./data','1976-2016_5+','cat',None,N=True)
    get_GCNN_data('./data','1976-2016_5+','solv',['cat'],N=True)
    get_GCNN_data('./data','1976-2016_5+','reag0',['cat','solv'],N=True)


    
