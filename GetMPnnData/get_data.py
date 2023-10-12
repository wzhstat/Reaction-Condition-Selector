import csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import csv
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

def remove_reagent(Smarts):
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
            

def save_MPNN_data(target,target_list):
    data = pd.read_csv('data/1976-2016_5+.csv')
    rxnsmile = Parallel(n_jobs=-1, verbose=4)(delayed(remove_reagent)(reaction) for reaction in list(data['reaction']))
    MLP_all_data = []
    l = len(rxnsmile)
    for i in range(l):
        if i % (l//100) == 0:
            print(i/(l//100),'%')
        dic = {}
        if data[target][i] not in target_list:
            continue
        dic['smarts'] = rxnsmile[i]
        dic[target] = target_list.index(data[target][i])
        MLP_all_data .append(dic)
    header = ['smarts',target]
    with open('data/GCN_%s_data_withoutN.csv'%target,'w',newline='') as f:
        f_csv = csv.DictWriter(f,header)
        f_csv.writeheader()
        f_csv.writerows(MLP_all_data)

def get_MPNN_data(path,withN,target_list):
    for target in  target_list:
        print('start to get %s data'%target)
        get_condition_list(target,path,withN)
        save_MPNN_data(target,target_list)
        print('get %s data done'%target)
        print('------------------')
   
