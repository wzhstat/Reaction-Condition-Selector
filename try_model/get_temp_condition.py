import csv
import pandas as pd
from rdkit import Chem
csv.field_size_limit(500 * 1024 * 1024)
def Cleam_merge(alist):
    newlist = eval(str(alist))
    outlist = []
    for r in newlist:
        mols = r.split('.')
        ol = []
        for mol in mols:
            try:
                mol = Chem.CanonSmiles(mol)
                ol.append(mol)
            except Exception as e:
                print(e)
        ol.sort()
        s = '.'.join(i for i in ol)
        if s not in outlist:
            outlist.append(s)
        outlist.sort()
    return outlist


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


all_data = []
data = pd.read_csv('data/data_5+.csv')
tem = data['template'][0]
adic = {}
adic['tem'] = tem
adic['cat'] = []
adic['solv'] = []
for i in range(len(data['template'])):
    if data['template'][i] != tem:
        all_data.append(adic)
        tem = data['template'][i]
        adic = {}
        adic['tem'] = tem
        adic['cat'] = []
        adic['solv'] = []
        if data['cat'][i] in cat_list and cat_list.index(data['cat'][i]) not in adic['cat']:
            adic['cat'].append(cat_list.index(data['cat'][i]))
        if data['solv'][i] in solv_list and solv_list.index(data['solv'][i]) not in adic['solv']:
            adic['solv'].append(solv_list.index(data['solv'][i]))
    else:
        if data['cat'][i] in cat_list and cat_list.index(data['cat'][i]) not in adic['cat']:
            adic['cat'].append(cat_list.index(data['cat'][i]))
        if data['solv'][i] in solv_list and solv_list.index(data['solv'][i]) not in adic['solv']:
            adic['solv'].append(solv_list.index(data['solv'][i]))

all_data = pd.DataFrame(all_data)
print(all_data)
print(len(all_data))
all_data.to_csv('./data/temp_condition.csv')
print('done')