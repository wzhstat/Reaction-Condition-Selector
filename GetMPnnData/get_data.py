import csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import csv
import sys
csv.field_size_limit(500 * 1024 * 1024)

def get_condition_list(target,path,withN = False):
    if withN:
        N = 'withN'
    else:
        N = 'withoutN'
    if target in ['cat','solv']:
        with open('%s/all_%s_%s.csv'%(path,target,withN),'r') as f:
            reader = csv.DictReader(f)
            for classes in reader:
                condition_list = list(classes.keys())
                print('%s num'%target,len(condition_list))
        return condition_list

    if target in ['reag0','reag1','reag2','reag3']:
        with open('%s/all_reag_%s.csv'%(path,withN),'r') as f:
            reader = csv.DictReader(f)
            for classes in reader:
                reag_list = list(classes.keys())
                print('reag num',len(reag_list))
        return reag_list
    
    else:
        raise KeyError("target must be in 'cat,solv,reag0,reag1,reag2,reag3'")

def remove_reagent(Smarts:str):
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
            

def open_MPNN_data(path: str,file_name: str,target: str,target_list: list, data: pd.DataFrame):
    data = data
    rxnsmile = Parallel(n_jobs=-1, verbose=4)(delayed(remove_reagent)(reaction) for reaction in list(data['reaction']))
    MLP_all_data = []
    l = len(rxnsmile)
    for i in range(l):
        if i % (l//100) == 0:
            print("\r", end="")
            print("Progress: {}%: ".format(i/(l//100)), "▋" * int(i/(l//50)), end="")
            sys.stdout.flush()
        dic = {}
        if data[target][i] not in target_list:
            continue
        dic['smarts'] = rxnsmile[i]
        dic[target] = target_list.index(data[target][i])
        MLP_all_data.append(dic)
    return MLP_all_data

def add_condition(data , condition:str, MLP_all_data: [dict]):
    l = len(MLP_all_data)
    for i in range(len(MLP_all_data)):
        if i % (l//100) == 0:
            print("\r", end="")
            print("Progress: {}%: ".format(i/(l//100)), "" * int(i/(l//50)), end="")
            sys.stdout.flush()
        MLP_all_data[i][condition] = data[condition][i]
    return MLP_all_data

def save_MPNN_csv(path,condition,MLP_all_data,target,N):
    if condition == None:
        header = ['smarts',target]
    else:
        header = ['smarts',condition,target]
    if N:
        path = '%s/GCN_%s_data_withN.csv'%(path,target)
    else:
        path = '%s/GCN_%s_data_withoutN.csv'%(path,target)
    with open(path,'w',newline='') as f:
        f_csv = csv.DictWriter(f,header)
        f_csv.writeheader()
        f_csv.writerows(MLP_all_data)

def get_GCNN_data(path: str, file_name: str, target: str, conditon: str, N: bool = False):
    data = pd.read_csv('%s/%s.csv'%(path,file_name))
    target_list = get_condition_list(target,path,withN = N)
    MLP_all_data = open_MPNN_data(path,file_name,target,target_list,data)
    if conditon != None:
        MLP_all_data = add_condition(data,conditon,MLP_all_data)
    save_MPNN_csv(path,conditon,MLP_all_data,target,N)

if __name__ == '__main__':
    get_GCNN_data('./data','1976-2016_5+','cat',None,N=True)
    get_GCNN_data('./data','1976-2016_5+','solv','cat',N=True)
    get_GCNN_data('./data','1976-2016_5+','reag0','solv',N=True)
    get_GCNN_data('./data','1976-2016_5+','cat',None)
    get_GCNN_data('./data','1976-2016_5+','solv','cat')
    get_GCNN_data('./data','1976-2016_5+','reag0','solv')
   
