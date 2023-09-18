import csv
import pandas as pd
from rdkit import Chem
csv.field_size_limit(500 * 1024 * 1024)
data = pd.read_csv('data/data.csv')

data_5 = []
n = 0
tem_num = 0
tem = data['tem_smart'][0]
tem_data = []
for i in range(len(data['tem_smart'])):
    adic = {}
    adic['reaction'] = data['reaction'][i]
    adic['products'] = data['products'][i]
    adic['reactants'] = data['reactants'][i]
    adic['tem_smart'] = data['tem_smart'][i]
    adic['template'] = n
    adic['cat'] = data['cat'][i]
    adic['solv'] = data['solv'][i]
    if tem != data['tem_smart'][i]:
        if tem_num >= 5 and tem != 'else':
            for j in tem_data:
                data_5.append(j)
            n +=1
        tem_num = 0
        tem_data = []
        tem = data['tem_smart'][i]
        print(n)
    else:
        tem_data.append(adic)
        tem_num += 1

all_data = pd.DataFrame(data_5)
print(all_data)
print(len(all_data))
all_data.to_csv('./data/data_5+.csv')
print('done')
    