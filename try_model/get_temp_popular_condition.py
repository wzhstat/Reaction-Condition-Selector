import csv
import pandas as pd
from rdkit import Chem
csv.field_size_limit(500 * 1024 * 1024)

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
bdic_cat = {}
bdic_solv = {}
adic['tem'] = tem
adic['cat'] = []
adic['solv'] = []
for i in range(len(data['template'])):
    if data['template'][i] != tem:
        bdic_cat = {k: v for k, v in sorted(bdic_cat.items(), key=lambda x: x[1], reverse=True)}
        bdic_solv = {k: v for k, v in sorted(bdic_solv.items(), key=lambda x: x[1], reverse=True)}
        cats = list(bdic_cat.keys())
        solvs = list(bdic_solv.keys())
        if len(cats) >= 10:
            adic['cat'] = cats[:10]
        else:
            adic['cat'] = cats+ ['None']*(10-len(cats))
        if len(solvs) >=10:
            adic['solv'] = solvs[:10]
        else:
            adic['solv'] = solvs + ['None']*(10-len(solvs))
        bdic_cat = {}
        bdic_solv = {}
        all_data.append(adic)
        tem = data['template'][i]
        adic = {}
        adic['tem'] = tem
        adic['cat'] = []
        adic['solv'] = []
        if data['cat'][i] in cat_list:
            if cat_list.index(data['cat'][i]) not in bdic_cat:
                bdic_cat[cat_list.index(data['cat'][i])] = 1
            else:
                bdic_cat[cat_list.index(data['cat'][i])] += 1
        if data['solv'][i] in solv_list:
            if solv_list.index(data['solv'][i]) not in bdic_solv:
                bdic_solv[solv_list.index(data['solv'][i])] = 1
            else:
                bdic_solv[solv_list.index(data['solv'][i])] += 1
    else:
        if data['cat'][i] in cat_list:
            if cat_list.index(data['cat'][i]) not in bdic_cat:
                bdic_cat[cat_list.index(data['cat'][i])] = 1
            else:
                bdic_cat[cat_list.index(data['cat'][i])] += 1
        if data['solv'][i] in solv_list:
            if solv_list.index(data['solv'][i]) not in bdic_solv:
                bdic_solv[solv_list.index(data['solv'][i])] = 1
            else:
                bdic_solv[solv_list.index(data['solv'][i])] += 1

all_data = pd.DataFrame(all_data)
print(all_data)
print(len(all_data))
all_data.to_csv('./data/temp_popular_condition.csv')
print('done')
