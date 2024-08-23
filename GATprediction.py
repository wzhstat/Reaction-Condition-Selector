from GAT.GAT_models import GAT
from GAT.get_dataset import GraphDataset, collate_fn, GraphDataLoader
import pandas as pd
import csv
import json,gzip
import argparse
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from joblib import Parallel, delayed
import torch
from tqdm import tqdm
import time
import os

def make_predictions(test_path,num_classes,model_weights_path, target = None):
    if torch.cuda.is_available():
        torch.cuda.set_device(3)
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')
    data_test = pd.read_csv(test_path)
    node_in_channels = 165
    edge_in_channels = 193
    hidden_channels = 1024
    num_classes = num_classes
    num_layers = 5
    num_heads = 8
    dropout =0.2
    model = GAT(node_in_channels, num_classes, num_layers, num_heads, dropout, edge_in_channels)
    model.load_state_dict(torch.load("%s.pth"%model_weights_path))
    test_dataset = GraphDataset(data_test,target)
    test_loader = GraphDataLoader(test_dataset, batch_size=1, collate_fn=collate_fn)
    pred = []
    model.eval()
    model.to(device)
    with torch.no_grad():
        for batch in tqdm(test_loader, desc="Processing batches"):
            batch = batch.to(device)
            outputs = model(batch)
            outputs = torch.softmax(outputs, dim=1)
            outputs.to(torch.device('cpu'))
            pred.append(outputs)
            torch.cuda.empty_cache()
    return pred
            

def encode_condition(condition_list:list,cat_list:list,solv_list:list,reag_list:list):
    '''
    encode condition list to index
    args:
        condition_list: condition list
        cat_list: catalyst list
        solv_list: solvent list
        reag_list: reagent list
    '''
    out = []
    for condition in condition_list:
        condition = eval(condition)
        cat = cat_list.index(condition[0])
        solv1 = solv_list.index(condition[1])
        solv2 = solv_list.index(condition[2])
        reag1 = reag_list.index(condition[3])
        reag2 = reag_list.index(condition[4])
        reag3 = reag_list.index(condition[5])
        out.append([cat,solv1,solv2,reag1,reag2,reag3])
    return out

def decode_condition(condition_list:list,cat_list:list,solv_list:list,reag_list:list):
    '''
    decode condition list to index
    args:
        condition_list: condition list
        cat_list: catalyst list
        solv_list: solvent list
        reag_list: reagent list
    '''
    out = []
    for condition in condition_list:
        cat = cat_list[condition[0]]
        solv1 = solv_list[condition[1]]
        solv2 = solv_list[condition[2]]
        reag1 = reag_list[condition[3]]
        reag2 = reag_list[condition[4]]
        reag3 = reag_list[condition[5]]
        out.append([cat,solv1,solv2,reag1,reag2,reag3])
    return out


def get_condition_key(path:str):
    with open('%s/cat_labels.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        cat_list_N = [row['cat'] for row in reader]

    with open('%s/solv_labels.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        solv_list_N = [row['solv'] for row in reader]

    with open('%s/reag_labels.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        reag_list_N = [row['reag'] for row in reader]   
    
    return cat_list_N,solv_list_N,reag_list_N


def get_condition_score(conditions,MPNN_out,condition_key):
    cat_list,solv_list,reag_list = condition_key
    condition_score = dict()
    for condition in conditions:
        text_condition = decode_condition([condition],cat_list,solv_list,reag_list)
        condition_score[str(text_condition[0])] = cal_condition_score(condition,MPNN_out)
    condition_score = sorted(condition_score.items(), key=lambda x:x[1],reverse=True)
    return condition_score

def cal_condition_score(condition,MPNN_out):
    score = 1
    for i in range(len(condition)-1):
        try:
            score *= MPNN_out[i][condition[i]]
        except:
            score *= 1e-10
    return float(score)

def condition_selector(args,template,MPNN_out, condition_library):
    try:
        cat_list,solv_list,reag_list = condition_key
        conditions = list(condition_library[template]['conditions'].keys())
        conditions = encode_condition(conditions,cat_list,solv_list,reag_list)
        condition_score = get_condition_score(conditions,MPNN_out,condition_key)
        return condition_score
    except:
        return []


def Prediction(args):
    '''
    This function is used to predict reaction conditions based on MPNN model,this function will give non-clustered results.
    args:
        args.test_path: path to test data
        args.model_path: path to model
        args.key_path: path to condition keys
        args.library_path: path to classed conditions library
        args.save_path: path to save condition prediction
    '''
    global condition_key
    t1 = time.time()
    # Load data
    test_data = pd.read_csv(args.test_path)
    template_r0 = test_data['tpl_SMARTS_r0']
    template_r1 = test_data['tpl_SMARTS_r1']
    template_r_1 = test_data['tpl_SMARTS_r0*']
    ids = test_data['_id']
    # MPNN prediction
    MPNN_pred = {}
    num_classes_dic = {"cat":439, "solv0":542, "solv1":542, "reag0":2746, "reag1":2746, "reag2":2746 }
    for target in ['cat','solv0','solv1','reag0','reag1','reag2']:
        model_dir = "%s/model_%s"%(args.model_path,target)
        MPNN_pred[target] = make_predictions(args.test_path,num_classes_dic[target],model_dir)
        print("done %s"%target)
    t2 = time.time()
    print('time:',t2-t1)
    # Load condition key
    condition_key = get_condition_key(args.label_path)

    # Load condition_library
    with gzip.open(args.library_path+'/condition_library_r0_1.json.gz','r') as f:
        conditions_library_r_1 = json.load(f)
    with gzip.open(args.library_path+'/condition_library_r0.json.gz','r') as f:
        conditions_library_r0 = json.load(f)
    with gzip.open(args.library_path+'/condition_library_r1.json.gz','r') as f:
        conditions_library_r1 = json.load(f)
    
    # Get condition prediction
    condition_pred = {}
    for i in range(test_data.shape[0]):
        if template_r1[i] in conditions_library_r1:
            condition_pred[str(ids[i])] = condition_selector(args,template_r1[i],[list(MPNN_pred['cat'][i])[0],list(MPNN_pred['solv0'][i])[0],list(MPNN_pred['solv1'][i])[0],list(MPNN_pred['reag0'][i])[0],list(MPNN_pred['reag1'][i])[0],list(MPNN_pred['reag2'][i])[0]],conditions_library_r1)
        elif template_r0[i] in conditions_library_r0:
            condition_pred[str(ids[i])] = condition_selector(args,template_r0[i],[list(MPNN_pred['cat'][i])[0],list(MPNN_pred['solv0'][i])[0],list(MPNN_pred['solv1'][i])[0],list(MPNN_pred['reag0'][i])[0],list(MPNN_pred['reag1'][i])[0],list(MPNN_pred['reag2'][i])[0]],conditions_library_r0)      
        else:
            condition_pred[str(ids[i])] = condition_selector(args,template_r_1[i],[list(MPNN_pred['cat'][i])[0],list(MPNN_pred['solv0'][i])[0],list(MPNN_pred['solv1'][i])[0],list(MPNN_pred['reag0'][i])[0],list(MPNN_pred['reag1'][i])[0],list(MPNN_pred['reag2'][i])[0]],conditions_library_r_1)
    t3 = time.time()
    print('time:',t3-t2)
    # Save
    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)
    with open('%s/condition_prediction.json'%args.save_path,'w') as f:
        json.dump(condition_pred,f)
    t4 = time.time()
    print('Save to: %s'%args.save_path)
    print(t4-t3)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='reaction condition prediction')
    parser.add_argument('--test_path', type=str, default='./data/data_test.csv', help='path to test data')
    parser.add_argument('--model_path', type=str, default='./GATmodel', help='path to model')
    parser.add_argument('--label_path', type=str, default='./data/labels', help='path to condition keys')
    parser.add_argument('--library_path', type=str, default='./data/condition_library', help='path to classed conditions library')
    parser.add_argument('--save_path', type=str, default='./data/GATprediction', help='path to save condition prediction')
    args = parser.parse_args()
    Prediction(args)
