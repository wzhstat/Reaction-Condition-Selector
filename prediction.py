from chemprop.args import PredictArgs
from chemprop.train.make_predictions import make_predictions
import pandas as pd
import csv
import json,gzip
import argparse
import time
import os
import torch


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


def get_condition_labels(path:str):
    '''
    get condition labels
    This function is used to get condition labels from csv files
    '''
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
    '''
    This function is used to obtain the scores of all candidates
    '''
    cat_list,solv_list,reag_list = condition_key
    condition_score = dict()
    for condition in conditions:
        text_condition = decode_condition([condition],cat_list,solv_list,reag_list)
        condition_score[str(text_condition[0])] = cal_condition_score(condition,MPNN_out)
    condition_score = sorted(condition_score.items(), key=lambda x:x[1],reverse=True)
    return condition_score

def cal_condition_score(condition,MPNN_out):
    '''
    The score of the complete reaction condition is the product of the probabilities of each of the reaction components
    '''
    score = 1
    for i in range(len(condition)-1):
        try:
            score *= MPNN_out[i][condition[i]]
        except:
            score *= 1e-10
    return score

def condition_selector(args,template,MPNN_out, condition_library):
    try:
        cat_list,solv_list,reag_list = condition_key
        conditions = list(condition_library[template]['conditions'].keys())
        conditions = encode_condition(conditions,cat_list,solv_list,reag_list)
        condition_score = get_condition_score(conditions,MPNN_out,condition_key)
        return condition_score
    except:
        return []

def MPNN_prediction(args,model_dir,smiles):
    MPNN_args = ['--test_path', '%s'%args.test_path, '--checkpoint_dir', '%s'%model_dir, '--preds_path', './sample_preds.csv']
    MPNN_args = PredictArgs().parse_args(MPNN_args)
    preds = make_predictions(MPNN_args,smiles)
    return preds

def Nu_condition_selector(MPNN_out,n_list):
    '''
    We use this naive candidate generation when the template for the predicted reaction is not in our template-condition library.
    The model captures the complete conditions directly by combining Top-3 predictions for each component
    '''
    cat_list,solv_list,reag_list = condition_key
    top3_indices = []
    for i in range(len(MPNN_out)):
        lists = torch.tensor(MPNN_out[i])
        top3_indices.append(torch.topk(lists, 3).indices)
    #combine the top-3 indices
    output = {}
    for i in range(n_list[0]):
        for j in range(n_list[1]):
            for k in range(n_list[2]):
                for l in range(n_list[3]):
                    for m in range(n_list[4]):
                        for n in range(1):
                            indices = [top3_indices[0][i], top3_indices[1][j], top3_indices[2][k], top3_indices[3][l], top3_indices[4][m], top3_indices[5][n]]
                            indices = [int(index) for index in indices]
                            score = 1
                            text_condition = decode_condition([indices],cat_list,solv_list,reag_list)
                            for con in range(len(indices)):
                                if score > 0:
                                    score *= MPNN_out[con][indices[con]]
                                else:
                                    score *= 1e-10
                            output[str(text_condition[0])] = score
    output = sorted(output.items(), key=lambda x:x[1],reverse=True)
    return output


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
    smiles = [[test_data['reaction'][i]] for i in range(len(test_data))]
    template_r0 = test_data['tpl_SMARTS_r0']
    template_r1 = test_data['tpl_SMARTS_r1']
    template_r_1 = test_data['tpl_SMARTS_r0*']
    ids = test_data['_id']
    # MPNN prediction
    MPNN_pred = {}
    for target in ['cat','solv0','solv1','reag0','reag1','reag2']:
        model_dir = "%s/MPNN_%s"%(args.model_path,target)
        MPNN_pred[target] = MPNN_prediction(args,model_dir,smiles)
    t2 = time.time()
    print('time:',t2-t1)
    # Load condition key
    condition_key = get_condition_labels(args.label_path)
    cat_list,solv_list,reag_list = condition_key

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
            condition_pred[str(ids[i])] = condition_selector(args,template_r1[i],[list(MPNN_pred['cat'][i][0]),list(MPNN_pred['solv0'][i][0]),list(MPNN_pred['solv1'][i][0]),list(MPNN_pred['reag0'][i][0]),list(MPNN_pred['reag1'][i][0]),list(MPNN_pred['reag2'][i][0])],conditions_library_r1)
        elif template_r0[i] in conditions_library_r0:
            condition_pred[str(ids[i])] = condition_selector(args,template_r0[i],[list(MPNN_pred['cat'][i][0]),list(MPNN_pred['solv0'][i][0]),list(MPNN_pred['solv1'][i][0]),list(MPNN_pred['reag0'][i][0]),list(MPNN_pred['reag1'][i][0]),list(MPNN_pred['reag2'][i][0])],conditions_library_r0)      
        elif template_r_1[i] in conditions_library_r_1:
            condition_pred[str(ids[i])] = condition_selector(args,template_r_1[i],[list(MPNN_pred['cat'][i][0]),list(MPNN_pred['solv0'][i][0]),list(MPNN_pred['solv1'][i][0]),list(MPNN_pred['reag0'][i][0]),list(MPNN_pred['reag1'][i][0]),list(MPNN_pred['reag2'][i][0])],conditions_library_r_1)
        else:
            condition_pred[str(ids[i])] = Nu_condition_selector([list(MPNN_pred['cat'][i][0]),list(MPNN_pred['solv0'][i][0]),list(MPNN_pred['solv1'][i][0]),list(MPNN_pred['reag0'][i][0]),list(MPNN_pred['reag1'][i][0]),list(MPNN_pred['reag2'][i][0])],[3,3,3,3,3])
    t3 = time.time()
    print('time:',t2-t1)
    # Save
    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)
    with open('%s/condition_prediction.json'%args.save_path,'w') as f:
        json.dump(condition_pred,f)
    t2 = time.time()
    print('Save to: %s'%args.save_path)
    print(t2-t1)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='reaction condition prediction')
    parser.add_argument('--test_path', type=str, default='./data/data_test.csv', help='path to test data')
    parser.add_argument('--model_path', type=str, default='./models', help='path to model')
    parser.add_argument('--label_path', type=str, default='./data/labels', help='path to condition keys')
    parser.add_argument('--library_path', type=str, default='./data/condition_library', help='path to classed conditions library')
    parser.add_argument('--save_path', type=str, default='./data/prediction', help='path to save condition prediction')
    args = parser.parse_args()
    Prediction(args)
