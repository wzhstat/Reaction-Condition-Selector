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

'''
with open('./data/all_cat3.csv','r') as f:
    reader = csv.DictReader(f)
    for classes in reader:
        cat_list = list(classes.keys())[:1183]
        print(len(cat_list))

with open('./data/all_solv3.csv','r') as f:
    reader = csv.DictReader(f)
    for classes in reader:
        solv_list = list(classes.keys())[:3646]
        print(len(solv_list))


cat_list.append('else')
solv_list.append('else')
'''
n = 0
all_data = []
data = pd.read_csv('data/1976-2016_3.csv')
with open('data/classif_by_temp2.csv','r') as f:
    reader = csv.DictReader(f)
    for classes in reader:
        for tem in classes:
            if n%338 == 0:
                print('已完成:',n/338,'%')
            go = False
            for i in eval(classes[tem]):
                adic = {}
                cats = data['catalyst'][int(i)-1]
                solvs = data['solvent'][int(i)-1]
                if pd.isna(cats):
                    cat = 'None'
                else:
                    cats = Cleam_merge(cats)
                    cat = '.'.join(i for i in cats)
                if pd.isna(solvs):
                    solv = 'None'
                else:
                    solvs = Cleam_merge(solvs)
                    solv = '.'.join(i for i in solvs)
                if cat == 'None' and solv == 'None':
                    continue
                adic['reaction'] = data['reaction'][int(i)-1]
                adic['products'] = '.'.join(i for i in eval(data['products'][int(i)-1]))
                adic['reactants'] = '.'.join(i for i in eval(data['reactants'][int(i)-1]))
                adic['template'] = n
                adic['tem_smart'] = tem
                adic['cat'] = cat
                adic['solv'] = solv
                go = True
                all_data .append(adic)
            if go:
                n +=1
                
all_data = pd.DataFrame(all_data)
print(all_data)
print(len(all_data))
all_data.to_csv('./data/data.csv')
print('done')

