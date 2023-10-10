import csv
import pandas as pd
from rdkit import Chem
csv.field_size_limit(500 * 1024 * 1024)


def Cleam_merge(alist):
    '''
    Clean and merge a list of dictionaries.

    Args:
        alist (list): A list of dictionaries.

    Returns:
        A merged list.
    '''
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

def count_data(data_name,save_path = "./data",m1 = 5,mc = 5,ms = 5, mr = 5):
    file_path= "%s/%s_%s+.csv"%(save_path,data_name,m1)
    '''
    Count the common data in a dataset.

    Args:
        data_name (str): The name of the dataset.
        save_path (str): The path to save the results.
        m1 (int): The minimum number of occurrences for a reaction type.
        mc (int): The minimum number of occurrences for a condition.
        ms (int): The minimum number of occurrences for a solvent.
        mr (int): The minimum number of occurrences for a reagent.

    '''
    n = 0
    # Load the data.
    data = pd.read_csv(file_path)
    tem = data['template'][0]
    all_cat_withN = {'None':0}
    all_cat_withoutN = {}
    all_solv_withN = {'None':0}
    all_solv_withoutN = {}
    all_reag_withN = {'None':0}
    all_reag_withoutN = {}
    lens = len(data['template'])

    # Count the common data.
    for i in range(lens):
        # Print the progress.
        if n% (lens//100) == 0:
            print(n/(lens//100),'%')
        
        # Get the data.
        cat = data['cat'][i]
        solv = data['solv'][i]
        reags = [data['reag%s'%j][i] for j in range(4)]

        # Count the common data.
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

    # Sort the results.
    all_cat_withN = {k: v for k, v in sorted(all_cat_withN.items(), key=lambda x: x[1], reverse=True)}
    all_cat_withoutN = {k: v for k, v in sorted(all_cat_withoutN.items(), key=lambda x: x[1], reverse=True)}
    all_solv_withN = {k: v for k, v in sorted(all_solv_withN.items(), key=lambda x: x[1], reverse=True)}
    all_solv_withoutN = {k: v for k, v in sorted(all_solv_withoutN.items(), key=lambda x: x[1], reverse=True)}
    all_reag_withN = {k: v for k, v in sorted(all_reag_withN.items(), key=lambda x: x[1], reverse=True)}
    all_reag_withoutN = {k: v for k, v in sorted(all_reag_withoutN.items(), key=lambda x: x[1], reverse=True)}
    all_cat_withN = {k: v for k, v in all_cat_withN.items() if v >= mc}
    all_cat_withoutN = {k: v for k, v in all_cat_withoutN.items() if v >= mc}
    all_solv_withN = {k: v for k, v in all_solv_withN.items() if v >= ms}
    all_solv_withoutN = {k: v for k, v in all_solv_withoutN.items() if v >= ms}
    all_reag_withN = {k: v for k, v in all_reag_withN.items() if v >= mr}
    all_reag_withoutN = {k: v for k, v in all_reag_withoutN.items() if v >= mr}


    # Save the results.    
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

