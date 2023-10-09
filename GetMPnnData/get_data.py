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

with open('./data/all_cat_withoutN.csv','r') as f:
    reader = csv.DictReader(f)
    for classes in reader:
        cat_list = list(classes.keys())
        print(len(cat_list))

with open('./data/all_solv_withoutN.csv','r') as f:
    reader = csv.DictReader(f)
    for classes in reader:
        solv_list = list(classes.keys())
        print(len(solv_list))

with open('./data/all_reag_withoutN.csv','r') as f:
    reader = csv.DictReader(f)
    for classes in reader:
        reag_list = list(classes.keys())
        print(len(reag_list))

with open('./data/all_reag_withN.csv','r') as f:
    reader = csv.DictReader(f)
    for classes in reader:
        reag_list1 = list(classes.keys())
        print(len(reag_list1))

def remove_reagent(Smarts):
    rxn = Smarts.split('>')
    reactant = rxn[0]
    product = rxn[2]
    reactant = reactant.split('.')
    product = product.split('.')
    outreactant = []
    for i in reactant:
        if 'C:' in i or 'CH:' in i or 'CH2:' in i or 'CH3:' in i:
            outreactant.append(i)
    rout = '.'.join(outreactant)
    outproduct = []
    for i in product:
        if 'C:' in i or 'CH:' in i or 'CH2:' in i or 'CH3:' in i:
            outproduct.append(i)
    pout = '.'.join(outproduct)
    out = rout + '>>' + pout
    return out
            

    


def get_MPNN_data(target,target_list):
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

data = pd.read_csv('data/1976-2016_5+.csv')
rxnsmile = Parallel(n_jobs=-1, verbose=4)(delayed(remove_reagent)(reaction) for reaction in list(data['reaction']))
get_MPNN_data('reag1',reag_list)