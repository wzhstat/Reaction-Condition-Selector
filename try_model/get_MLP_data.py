from rdkit.Chem import AllChem
import rdkit.Chem as Chem
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import csv

def get_rxnfp(reaction):
    (reactant,product) = reaction
    rm = Chem.MolFromSmiles(reactant)
    pm = Chem.MolFromSmiles(product)
    info = {}
    rfpgen= np.array(AllChem.GetMorganFingerprintAsBitVect(rm, useChirality=True, radius=2, nBits = 512, bitInfo=info))
    pfpgen= np.array(AllChem.GetMorganFingerprintAsBitVect(pm, useChirality=True, radius=2, nBits = 512, bitInfo=info))
    rxnfp = pfpgen-rfpgen
    return (rfpgen,pfpgen,rxnfp)

def get_b_tem(template):
    num = int(template)
    b_num = []
    while num != 0:
        b_num.append(num%2)
        num = num//2
    add_0 = [0]*(16-len(b_num))
    return add_0 + b_num

with open('./data/all_cat3.csv','r') as f:
    reader = csv.DictReader(f)
    for classes in reader:
        cat_list = list(classes.keys())[:866]
        print(len(cat_list))

with open('./data/all_solv3.csv','r') as f:
    reader = csv.DictReader(f)
    for classes in reader:
        solv_list = list(classes.keys())[:2702]
        print(len(solv_list))

'''

data = pd.read_csv('data/data_5+.csv')
rxnfps = Parallel(n_jobs=-1, verbose=4)(delayed(get_rxnfp)(reaction) for reaction in list(data[['reactants','products']].apply(tuple, axis=1)))
#b_num = Parallel(n_jobs=-1, verbose=4)(delayed(get_b_tem)(tem_class) for tem_class in list(data['template']))
#one_hot_temp =np.array(pd.get_dummies(data['template'],columns=data.columns))
MLP_cat_data = []
for i in range(len(rxnfps)):
    if i % 16890 == 0:
        print(i/16890,'%')
    dic = {}
    if data['cat'][i] not in cat_list:
        continue
    dic['rfpgen'] = rxnfps[i][0]
    dic['pfpgen'] = rxnfps[i][1]
    dic['rxnfp'] = rxnfps[i][2]
    dic['tem'] = data['template'][i] 
    dic['cat'] = cat_list.index(data['cat'][i])
    
    MLP_cat_data.append(dic)


MLP_cat_data = pd.DataFrame(MLP_cat_data)
print(MLP_cat_data)
MLP_cat_data.to_json('data/MLP_cat_data.json.gz', orient='records', compression='gzip')
print('done_cat')

MLP_solv_data = []
for i in range(len(rxnfps)):
    if i % 16890 == 0:
        print(i/16890,'%')
    dic = {}
    if data['solv'][i] not in solv_list:
        continue
    dic['rfpgen'] = rxnfps[i][0]
    dic['pfpgen'] = rxnfps[i][1]
    dic['rxnfp'] = rxnfps[i][2]
    dic['tem'] = data['template'][i] 
    dic['solv'] =solv_list.index(data['solv'][i])
    MLP_solv_data.append(dic)

MLP_solv_data = pd.DataFrame(MLP_solv_data)
print(MLP_solv_data)
MLP_solv_data.to_json('data/MLP_solv_data.json.gz', orient='records', compression='gzip')
print('done_solv')
'''

cat_list.append('None')
solv_list.append('None')
cat_list.append('else')
solv_list.append('else')
data = pd.read_csv('data/data_5+.csv')
MLP_all_data = []
rxnfps = Parallel(n_jobs=-1, verbose=4)(delayed(get_rxnfp)(reaction) for reaction in list(data[['reactants','products']].apply(tuple, axis=1)))
for i in range(len(rxnfps)):
    if i % 9485 == 0:
        print(i/9485,'%')
    dic = {}
    dic['rfpgen'] = rxnfps[i][0]
    dic['pfpgen'] = rxnfps[i][1]
    dic['rxnfp'] = rxnfps[i][2]
    dic['tem'] = data['template'][i] 
    if data['cat'][i] not in cat_list:
        dic['cat'] = cat_list.index('else')
    else:
        dic['cat'] = cat_list.index(data['cat'][i])
    if data['solv'][i] not in solv_list:
        dic['solv'] = solv_list.index('else')
    else:
        dic['solv'] =solv_list.index(data['solv'][i])
    MLP_all_data .append(dic)


MLP_all_data = pd.DataFrame(MLP_all_data)
MLP_all_data.to_json('data/MLP_all_data.json.gz', orient='records', compression='gzip')