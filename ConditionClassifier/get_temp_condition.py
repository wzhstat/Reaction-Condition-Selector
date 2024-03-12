import pandas as pd
import numpy as np
import json
import sys
import csv
from joblib import Parallel, delayed
import gzip
from func_timeout import func_timeout, FunctionTimedOut
from rdchiral import template_extractor

def get_condition_key(path:str):
    with open('%s/all_cat_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        cat_list_N = list(row['cat'] for row in reader)

    with open('%s/all_solv_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        solv_list_N = [row['solv'] for row in reader]

    with open('%s/all_reag_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        reag_list_N = [row['reag'] for row in reader]
    return cat_list_N,solv_list_N,reag_list_N


def get_temp_condition(args):
    cat_list_N,solv_list_N,reag_list_N = get_condition_key(args.save_path)
    print('start to get condition list')
    if args.data_set == 'all':
        data =pd.read_csv('%s/%s_%s+.csv'%(args.save_path,args.data_name,args.min_num_covered_rxns_by_rxn_centralized_template))
    elif args.data_set == 'train':
        data =pd.read_csv('%s/data_train.csv'%(args.save_path))
    elif args.data_set == 'test':
        data =pd.read_csv('%s/data_test.csv'%(args.save_path))
    elif args.data_set == 'val':
        data =pd.read_csv('%s/data_val.csv'%(args.save_path))
    data = data.sort_values(by=['template_%s'%args.tpl_radius])
    l = data.shape[0]
    condition_list = []
    adic = {'tpl':data['template_0'][0],'tpl_smarts':data['tpl_smarts_%s'%args.tpl_radius][0],'conditions':{}}
    for i in range(l):
        if i%(l//1000) == 0:
            print("\r", end="")
            print("Progress: {}%: ".format(i/((l//1000)*10)), "▋" * int(i/((l//1000)*20)), end="")
            sys.stdout.flush()
        tem = data['template_0'][i]
        cat = data['cat'][i] if str(data['cat'][i]) != 'nan' else 'None'
        solv0 = data['solv0'][i] if str(data['solv0'][i]) != 'nan' else 'None'
        solv1 = data['solv1'][i] if str(data['solv1'][i]) != 'nan' else 'None'
        reag0 = data['reag0'][i] if str(data['reag0'][i]) != 'nan' else 'None'
        reag1 = data['reag1'][i] if str(data['reag1'][i]) != 'nan' else 'None'
        reag2 = data['reag2'][i] if str(data['reag2'][i]) != 'nan' else 'None'
        if tem == adic['tpl']:
            if cat not in cat_list_N:
                continue
            elif solv1 not in solv_list_N or solv0 not in solv_list_N:
                continue
            elif reag2 not in reag_list_N or reag1 not in reag_list_N or reag0 not in reag_list_N:
                continue
            condition = [cat,solv0,solv1,reag0,reag1,reag2]
            for j in range(len(condition)):
                if str(condition[j]) == 'nan':
                    condition[j] = 'None'
            if str(condition) in adic['conditions']:
                adic['conditions'][str(condition)] += 1
            else:
                adic['conditions'][str(condition)] = 1
        else:
            condition_list.append(adic)
            for k in range(adic['tpl']+1,tem):
                condition_list.append({'tpl':k,'tpl_smarts':'None','conditions':{}})
            adic = {'tpl':tem,'tpl_smarts':data['tpl_smarts_%s'%args.tpl_radius][i],'conditions':{}}
            if cat not in cat_list_N:
                continue
            elif solv1 not in solv_list_N or solv0 not in solv_list_N:
                continue
            elif reag2 not in reag_list_N or reag1 not in reag_list_N or reag0 not in reag_list_N:
                continue
            condition = [cat,solv0,solv1,reag0,reag1,reag2]
            for j in range(len(condition)):
                if str(condition[j]) == 'nan':
                    condition[j] = 'None'
            if str(condition) in adic['conditions']:
                adic['conditions'][str(condition)] += 1
            else:
                adic['conditions'][str(condition)] = 1
    condition_list.append(adic)
    print('get condition list done')
    return condition_list

if __name__ == '__main__':
    #get_temp_condition()
    pass

