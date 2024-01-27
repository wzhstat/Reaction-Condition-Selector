import csv
import pandas as pd
from rdkit import Chem
import sys
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

def count_data(args):
    file_path= "%s/%s.csv"%(args.save_path,args.data_name)
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
    all_cat_withN = {'None':0}
    all_cat_withoutN = {}
    all_solv_withN = {'None':0}
    all_solv_withoutN = {}
    all_reag_withN = {'None':0}
    all_reag_withoutN = {}
    lens = data.shape[0]

    # Count the common data.
    for i in range(lens):
        # Print the progress.
        if n% (lens//100) == 0:
            print("\r", end="")
            print("Progress: {}%: ".format(n/(lens//100)), "▋" * int(0.5*n/(lens//100)), end="")
            sys.stdout.flush()
        
        # Get the data.
        cats = data['catalyst'][i]
        solvs = data['solvent'][i]
        reag = data['reagent'][i]
        if cats == 'None':
            cat = cats
        else:
            cats = Cleam_merge(cats)
            cat = '.'.join(i for i in cats)
        if solvs == 'None':
            solvs = [solvs]
        else:
            solvs = Cleam_merge(solvs)
        if reag == 'None':
            reags = [reag]
        else:
            reags = Cleam_merge(reag)

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
        for solv in solvs:
            if solv == 'None':
                all_solv_withN['None'] += 1
            else:
                if solv in all_solv_withoutN:
                    all_solv_withN[solv] += 1
                    all_solv_withoutN[solv] += 1
                else:
                    all_solv_withN[solv] = 1
                    all_solv_withoutN[solv] = 1
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
        n += 1

    # Sort the results.
    all_cat_withN = {k: v for k, v in sorted(all_cat_withN.items(), key=lambda x: x[1], reverse=True)}
    all_cat_withoutN = {k: v for k, v in sorted(all_cat_withoutN.items(), key=lambda x: x[1], reverse=True)}
    all_solv_withN = {k: v for k, v in sorted(all_solv_withN.items(), key=lambda x: x[1], reverse=True)}
    all_solv_withoutN = {k: v for k, v in sorted(all_solv_withoutN.items(), key=lambda x: x[1], reverse=True)}
    all_reag_withN = {k: v for k, v in sorted(all_reag_withN.items(), key=lambda x: x[1], reverse=True)}
    all_reag_withoutN = {k: v for k, v in sorted(all_reag_withoutN.items(), key=lambda x: x[1], reverse=True)}
    all_cat_withN = {k: v for k, v in all_cat_withN.items() if v >= args.min_num_covered_rxns_by_catalyst}
    all_cat_withoutN = {k: v for k, v in all_cat_withoutN.items() if v >= args.min_num_covered_rxns_by_catalyst}
    all_solv_withN = {k: v for k, v in all_solv_withN.items() if v >= args.min_num_covered_rxns_by_solvent}
    all_solv_withoutN = {k: v for k, v in all_solv_withoutN.items() if v >= args.min_num_covered_rxns_by_solvent}
    all_reag_withN = {k: v for k, v in all_reag_withN.items() if v >= args.min_num_covered_rxns_by_reagent}
    all_reag_withoutN = {k: v for k, v in all_reag_withoutN.items() if v >= args.min_num_covered_rxns_by_reagent}
    print(all_cat_withN)
    print(all_solv_withN)
    print(all_reag_withN)

    # Save the results.    
    with open('%s/all_cat_withN.csv'%args.save_path,'w', newline='') as catfile:
        header = ['cat','count']
        writer = csv.writer(catfile)
        writer.writerow(header)
        for key,value in all_cat_withN.items():
            writer.writerow([key,value])
       

    with open('%s/all_cat_withoutN.csv'%args.save_path,'w', newline='') as catfile:
        header = ['cat','count']
        writer = csv.writer(catfile)
        writer.writerow(header)
        for key,value in all_cat_withoutN.items():
            writer.writerow([key,value])

    with open('%s/all_solv_withN.csv'%args.save_path,'w', newline='') as solvfile:
        header = ['solv','count']
        writer = csv.writer(solvfile)
        writer.writerow(header)
        for key,value in all_solv_withN.items():
            writer.writerow([key,value])

    with open('%s/all_solv_withoutN.csv'%args.save_path,'w', newline='') as solvfile:
        header = ['solv','count']
        writer = csv.writer(solvfile)
        writer.writerow(header)
        for key,value in all_solv_withoutN.items():
            writer.writerow([key,value])
    
    with open('%s/all_reag_withN.csv'%args.save_path,'w', newline='') as reagfile:
        header = ['reag','count']
        writer = csv.writer(reagfile)
        writer.writerow(header)
        for key,value in all_reag_withN.items():
            writer.writerow([key,value])
    
    with open('%s/all_reag_withoutN.csv'%args.save_path,'w', newline='') as reagfile:
        header = ['reag','count']
        writer = csv.writer(reagfile)
        writer.writerow(header)
        for key,value in all_reag_withoutN.items():
            writer.writerow([key,value])

if __name__ == '__main__':
    pass
