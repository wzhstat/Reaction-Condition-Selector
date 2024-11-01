import pandas as pd
import json
from joblib import Parallel, delayed
import csv
import gzip
from operator import add
from functools import reduce
from ConditionClassifier import get_labeled_condition,same_class
import argparse
import numpy as np
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def remove_chirality(SMILES:str):
    '''
    This function is used to remove chirality from the reaction conditions
    '''
    SMILES = SMILES.replace('@','')
    SMILES = SMILES.replace('/','')
    SMILES = SMILES.replace('\\','')
    return SMILES

def get_true_data(path:str):
    data = pd.read_csv(path, dtype={'reag0':'str','reag1':'str','reag2':'str'})
    l = data.shape[0]
    y_true = []
    for i in range(l):
        condition_trun = []
        for i in [str(data['cat'][i]),str(data['solv0'][i]),str(data['solv1'][i]),str(data['reag0'][i]),str(data['reag1'][i]),str(data['reag2'][i])]:
            try:
                mol = Chem.MolFromSmiles(i)
                inch = Chem.MolToInchi(mol)
                mol = Chem.MolFromInchi(inch)
                i = Chem.MolToSmiles(mol)
            except:
                i = i
            if i == 'nan':
                i = 'None'
            condition_trun.append(i)
        y_true.append(condition_trun)
    return y_true

def get_condition_set(alist:list):
    return set('.'.join(alist).split('.'))
    

def is_topk_acc(y_true:list,pred:list):
    '''
    Determine if ground_truth is in the Top-K prediction.
    '''
    out = []
    y_true = [remove_chirality(x) for x in y_true]
    pred = [(remove_chirality(y) for y in x) for x in pred]
    y_true = get_condition_set(y_true)
    pred = [get_condition_set(x) for x in pred]
    for k in [1,3,5,10]:
        for i in pred[:k]:
            if i == y_true:
                out.append(1)
                break
        else:
            out.append(0)
    return out
        
    

def is_topk_cluster_acc(y_true:list,pred:list):
    '''
    Determine if ground_truth is in the Top-K prediction.
    '''
    out = []
    y_true_label = get_labeled_condition([y_true])[0][1]
    pred_label = [x[1] for x in get_labeled_condition(pred)]
    for k in [1,3,5,10]:
        for i in pred_label[:k]:
            if same_class(i,y_true_label):
                out.append(1)
                break
        else:
            out.append(0)
    return out


def cal_acc(args):
    '''
    This function calculates the accuracy of a model or algorithm. 
    args:
        args.data_path: str, the path to the test data.
        args.pred_path: str, the path to the prediction file.
    '''
    with open('%s/cluster_condition_prediction.json'%args.pred_path,'r') as f:
        cluster_pred = json.load(f)
    with open('%s/condition_prediction.json'%args.pred_path,'r') as f:
        pred = json.load(f)
    y_true = get_true_data(args.data_path)
    test_data_key = pd.read_csv(args.data_path, dtype={'reag0':'str','reag1':'str','reag2':'str'})['_id'].tolist()
    acc_list = Parallel(n_jobs=10)(delayed(is_topk_acc)(y_true[i],[eval(x[0]) for x in pred[str(test_data_key[i])][:10]]) for i in range(len(y_true)))
    acc_list = np.array(acc_list)
    clu_acc_list = Parallel(n_jobs=10)(delayed(is_topk_cluster_acc)(y_true[i],[x['best condition'] for x in cluster_pred[str(test_data_key[i])][1][:10]]) for i in range(len(y_true)))
    clu_acc_list = np.array(clu_acc_list)
    print('Top-1 Accuracy:',np.mean(acc_list[:,0]))
    print('Top-3 Accuracy:',np.mean(acc_list[:,1]))
    print('Top-5 Accuracy:',np.mean(acc_list[:,2]))
    print('Top-10 Accuracy:',np.mean(acc_list[:,3]))
    print('Cluster Top-1 Accuracy:',np.mean(clu_acc_list[:,0]))
    print('Cluster Top-3 Accuracy:',np.mean(clu_acc_list[:,1]))
    print('Cluster Top-5 Accuracy:',np.mean(clu_acc_list[:,2]))
    print('Cluster Top-10 Accuracy:',np.mean(clu_acc_list[:,3]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path',type=str,default='./data/data_test.csv')
    parser.add_argument('--pred_path',type=str,default='./data/prediction')
    args = parser.parse_args()
    cal_acc(args)

if __name__ == '__main__':
    main()

    



    


