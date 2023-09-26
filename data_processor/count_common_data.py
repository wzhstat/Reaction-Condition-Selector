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

def count_data(save_path = "./data",m = 5):
    file_path= "%s/data_%s+.csv"%(save_path,m)
    '''
    save_path: the floder path to save the data.
    This function is used to statistics the conditions.
    '''
    n = 0
    data = pd.read_csv(file_path)
    tem = data['template'][0]
    all_cat_withN = {'None':0}
    all_cat_withoutN = {}
    all_solv_withN = {'None':0}
    all_solv_withoutN = {}
    all_reag_withN = {'None':0}
    all_reag_withoutN = {}
    lens = len(data['template'])
    for i in range(lens):
        if n% (lens//100) == 0:
            print(n/(lens//100),'%')
        cat = data['cat'][i]
        solv = data['solv'][i]
        reags = [data['reag%s'%j][i] for j in range(4)]
        if cat == 'None':
            all_cat_withN['None'] += 1
        else:
            if cat in all_cat_withoutN:
                all_cat_withoutN[cat] += 1
                all_cat_withN[cat] += 1
            else:
                all_cat_withoutN[cat] = 1
                all_cat_withN[cat] = 1
        if solv == 'None':
            all_solv_withN['None'] += 1
        else:
            if solv in all_solv_withoutN:
                all_solv_withN[solv] += 1
                all_solv_withoutN[solv] += 1
            else:
                all_solv_withN[solv] = 1
                all_solv_withoutN[solv] = 1
        n += 1
        for reag in reags:
            if reag == 'None':
                all_reag_withN['None'] += 1
            else:
                if reag in all_reag_withoutN:
                    all_reag_withN[reag] += 1
                    all_reag_withoutN[reag] += 1
                else:
                    all_reag_withN[reag] = 1
                    all_reag_withoutN[reag] = 1

        
    all_cat_withN = {k: v for k, v in sorted(all_cat_withN.items(), key=lambda x: x[1], reverse=True)}
    all_cat_withoutN = {k: v for k, v in sorted(all_cat_withoutN.items(), key=lambda x: x[1], reverse=True)}
    all_solv_withN = {k: v for k, v in sorted(all_solv_withN.items(), key=lambda x: x[1], reverse=True)}
    all_solv_withoutN = {k: v for k, v in sorted(all_solv_withoutN.items(), key=lambda x: x[1], reverse=True)}
    all_reag_withN = {k: v for k, v in sorted(all_reag_withN.items(), key=lambda x: x[1], reverse=True)}
    all_reag_withoutN = {k: v for k, v in sorted(all_reag_withoutN.items(), key=lambda x: x[1], reverse=True)}


            
    with open('%s/all_cat_withN.csv'%save_path,'w', newline='') as catfile:
        writer = csv.writer(catfile)
        writer.writerow(all_cat_withN.keys())
        writer.writerow(all_cat_withN.values())

    with open('%s/all_cat_withoutN.csv'%save_path,'w', newline='') as catfile:
        writer = csv.writer(catfile)
        writer.writerow(all_cat_withoutN.keys())
        writer.writerow(all_cat_withoutN.values())

    with open('%s/all_solv_withN.csv'%save_path,'w', newline='') as solvfile:
        writer = csv.writer(solvfile)
        writer.writerow(all_solv_withN.keys())
        writer.writerow(all_solv_withN.values())

    with open('%s/all_solv_withoutN.csv'%save_path,'w', newline='') as solvfile:
        writer = csv.writer(solvfile)
        writer.writerow(all_solv_withoutN.keys())
        writer.writerow(all_solv_withoutN.values())
    
    with open('%s/all_reag_withN.csv'%save_path,'w', newline='') as reagfile:
        writer = csv.writer(reagfile)
        writer.writerow(all_reag_withN.keys())
        writer.writerow(all_reag_withN.values())
    
    with open('%s/all_reag_withoutN.csv'%save_path,'w', newline='') as reagfile:
        writer = csv.writer(reagfile)
        writer.writerow(all_reag_withoutN.keys())
        writer.writerow(all_reag_withoutN.values())


