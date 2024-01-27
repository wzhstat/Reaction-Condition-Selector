import pandas as pd
import numpy as np
import json
import sys
import csv
from joblib import Parallel, delayed
import gzip



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
    data.sort_values(by=['template'],inplace=True)
    l = data.shape[0]
    condition_list = []
    adic = {'tem':data['template'][0],'tem_smart':data['tem_smart'][0],'conditions':{}}
    for i in range(l):
        if i%(l//1000) == 0:
            print("\r", end="")
            print("Progress: {}%: ".format(i/((l//1000)*10)), "▋" * int(i/((l//1000)*20)), end="")
            sys.stdout.flush()
        tem = data['template'][i]
        cat = data['cat'][i] if str(data['cat'][i]) != 'nan' else 'None'
        solv0 = data['solv0'][i] if str(data['solv0'][i]) != 'nan' else 'None'
        solv1 = data['solv1'][i] if str(data['solv1'][i]) != 'nan' else 'None'
        reag0 = data['reag0'][i] if str(data['reag0'][i]) != 'nan' else 'None'
        reag1 = data['reag1'][i] if str(data['reag1'][i]) != 'nan' else 'None'
        reag2 = data['reag2'][i] if str(data['reag2'][i]) != 'nan' else 'None'
        if tem == adic['tem']:
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
            for k in range(adic['tem']+1,tem):
                condition_list.append({'tem':k,'tem_smart':'None','conditions':{}})
            adic = {'tem':tem,'tem_smart':data['tem_smart'][i],'conditions':{}}
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


    condition_list = pd.DataFrame(condition_list)
    condition_list.to_json('%s/condition_list_%s.json.gz'%(args.save_path,args.data_set),orient='records',compression='gzip')
    print('get condition list done')
    return condition_list

if __name__ == '__main__':
    #get_temp_condition()
    with gzip.open('data/condition_list.json.gz') as f:
        condition_list = json.load(f)
    print(condition_list[0])
    '''
    condition_list = classify_condiition_in_temp.get_tem_condition('data')
    classed_condition_list = list(Parallel(n_jobs=-1)(delayed(classify_condiition_in_temp.Classify_reaction_conditions)(sorted(condition_list['conditions'][i]),condition_list['tem'][i],condition_list['tem_smart'][i]) for i in range(len(condition_list))))
    with open('data/classed_condition_list.json','w') as f:
        json.dump(classed_condition_list,f)
    '''
