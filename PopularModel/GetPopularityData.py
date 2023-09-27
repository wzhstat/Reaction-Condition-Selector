import csv
import pandas as pd
from rdkit import Chem
import hashlib

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


def get_popularity_data(file_name = '1976-2016'):
    '''
    This function is used to get the most popular catalysts and solvents in an reaction template.
    '''
    n = 0
    all_data = []
    data = pd.read_csv('data/%s.csv'%(file_name))
    with open('data/classif_by_temp.csv','r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            for tem in classes:
                if tem == 'else':
                    continue
                if n%378 == 0:
                    print('已完成:',n/378,'%')
                n +=1
                catdic = {'None':0}
                solvdic = {'None':0}
                reagdic = {'None':0}
                for i in eval(classes[tem]):
                    cat = data['catalyst'][int(i)-1]
                    solv = data['solvent'][int(i)-1]
                    reag = data['reagent'][int(i)-1]
                    if cat == 'None':
                        catdic['None'] += 1
                    else:
                        cats = Cleam_merge(cat)
                        cat = '.'.join(i for i in cats)
                        if cat not in catdic:
                            catdic[cat] = 1
                        else:
                            catdic[cat] += 1
                    if solv == 'None':
                        solvdic['None'] += 1
                    else:
                        solvs = Cleam_merge(solv)
                        solv = '.'.join(i for i in solvs)
                        if solv not in solvdic:
                            solvdic[solv] = 1
                        else:
                            solvdic[solv] += 1
                    if reag == 'None':
                        reagdic['None'] += 1
                    else:
                        reags = eval(str(reag))
                        for reag in reags:
                            if reag not in reagdic:
                                reagdic[reag] = 1
                            else:
                                reagdic[reag] += 1  
                catdic = {k: v for k, v in sorted(catdic.items(), key=lambda x: x[1], reverse=True)}
                solvdic = {k: v for k, v in sorted(solvdic.items(), key=lambda x: x[1], reverse=True)}
                reagdic = {k: v for k, v in sorted(reagdic.items(), key=lambda x: x[1], reverse=True)}
                catkeys = list(catdic.keys())     #分别储存
                solvkeys = list(solvdic.keys())
                reagkeys = list(reagdic.keys())
                dic = {}
                dic['_id'] = n
                dic['Templates'] = tem
                if tem != 'else':
                    dic['hash_id'] = hashlib.md5(tem.encode()).hexdigest()
                else:
                    dic['hash_id'] = '0'
                m1 = 0
                while m1 <10:
                    name = 'cat_top-%d'%(m1+1)
                    if len(catdic) > m1:
                        dic[name] = catkeys[m1]
                    else:
                        dic[name] = 'None'
                    m1 += 1
                m2 = 0
                while m2 <10:
                    name = 'solv_top-%d'%(m2+1)
                    if len(solvdic) > m2:
                        dic[name] = solvkeys[m2]
                    else:
                        dic[name] = 'None'
                    m2 +=1
                m3 = 0
                while m3 <10:
                    name = 'reag_top-%d'%(m3+1)
                    if len(reagdic) > m3:
                        dic[name] = reagkeys[m3]
                    else:
                        dic[name] = 'None'
                    m3 +=1
                all_data.append(dic)
    all_data = sorted(all_data,key = lambda x: x['hash_id'])
    all_data = pd.DataFrame(all_data)
    all_data.to_csv('data/popularity_data.csv')
    print('done')


get_popularity_data()