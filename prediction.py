from chemprop.args import PredictArgs
from typing import List, Optional, Union, Tuple
from chemprop.train import load_model, set_features, load_data, predict_and_save
from chemprop.args import PredictArgs, TrainArgs
from chemprop.data import get_data,StandardScaler,MoleculeDataLoader, AtomBondScaler
from chemprop.models import MoleculeModel
from chemprop.uncertainty import UncertaintyCalibrator, build_uncertainty_calibrator
import pandas as pd
import csv
import json,gzip
import argparse
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from joblib import Parallel, delayed
import time
import os

def make_predictions(
    args: PredictArgs,
    smiles: List[List[str]] = None,
    model_objects: Tuple[
        PredictArgs,
        TrainArgs,
        List[MoleculeModel],
        List[Union[StandardScaler, AtomBondScaler]],
        int,
        List[str],
    ] = None,
    calibrator: UncertaintyCalibrator = None,
    return_invalid_smiles: bool = True,
    return_index_dict: bool = False,
    return_uncertainty: bool = False,
) -> List[List[Optional[float]]]:
    """
    Loads data and a trained model and uses the model to make predictions on the data.

    If SMILES are provided, then makes predictions on smiles.
    Otherwise makes predictions on :code:`args.test_data`.

    :param args: A :class:`~chemprop.args.PredictArgs` object containing arguments for
                loading data and a model and making predictions.
    :param smiles: List of list of SMILES to make predictions on.
    :param model_objects: Tuple of output of load_model function which can be called separately outside this function. Preloaded model objects should have
                used the non-generator option for load_model if the objects are to be used multiple times or are intended to be used for calibration as well.
    :param calibrator: A :class: `~chemprop.uncertainty.UncertaintyCalibrator` object, for use in calibrating uncertainty predictions.
                Can be preloaded and provided as a function input or constructed within the function from arguments. The models and scalers used
                to initiate the calibrator must be lists instead of generators if the same calibrator is to be used multiple times or
                if the same models and scalers objects are also part of the provided model_objects input.
    :param return_invalid_smiles: Whether to return predictions of "Invalid SMILES" for invalid SMILES, otherwise will skip them in returned predictions.
    :param return_index_dict: Whether to return the prediction results as a dictionary keyed from the initial data indexes.
    :param return_uncertainty: Whether to return uncertainty predictions alongside the model value predictions.
    :return: A list of lists of target predictions. If returning uncertainty, a tuple containing first prediction values then uncertainty estimates.
    """
    if model_objects:
        (
            args,
            train_args,
            models,
            scalers,
            num_tasks,
            task_names,
        ) = model_objects
    else:
        (
            args,
            train_args,
            models,
            scalers,
            num_tasks,
            task_names,
        ) = load_model(args, generator=True)

    num_models = len(args.checkpoint_paths)

    set_features(args, train_args)

    # Note: to get the invalid SMILES for your data, use the get_invalid_smiles_from_file or get_invalid_smiles_from_list functions from data/utils.py
    full_data, test_data, test_data_loader, full_to_valid_indices = load_data(
        args, smiles
    )

    if args.uncertainty_method is None and (args.calibration_method is not None or args.evaluation_methods is not None):
        if args.dataset_type in ['classification', 'multiclass']:
            args.uncertainty_method = 'classification'
        else:
            raise ValueError('Cannot calibrate or evaluate uncertainty without selection of an uncertainty method.')


    if calibrator is None and args.calibration_path is not None:

        calibration_data = get_data(
            path=args.calibration_path,
            smiles_columns=args.smiles_columns,
            target_columns=task_names,
            args=args,
            features_path=args.calibration_features_path,
            features_generator=args.features_generator,
            phase_features_path=args.calibration_phase_features_path,
            atom_descriptors_path=args.calibration_atom_descriptors_path,
            bond_descriptors_path=args.calibration_bond_descriptors_path,
            max_data_size=args.max_data_size,
            loss_function=args.loss_function,
        )

        calibration_data_loader = MoleculeDataLoader(
            dataset=calibration_data,
            batch_size=args.batch_size,
            num_workers=args.num_workers,
        )

        if isinstance(models, List) and isinstance(scalers, List):
            calibration_models = models
            calibration_scalers = scalers
        else:
            calibration_model_objects = load_model(args, generator=True)
            calibration_models = calibration_model_objects[2]
            calibration_scalers = calibration_model_objects[3]

        calibrator = build_uncertainty_calibrator(
            calibration_method=args.calibration_method,
            uncertainty_method=args.uncertainty_method,
            interval_percentile=args.calibration_interval_percentile,
            regression_calibrator_metric=args.regression_calibrator_metric,
            calibration_data=calibration_data,
            calibration_data_loader=calibration_data_loader,
            models=calibration_models,
            scalers=calibration_scalers,
            num_models=num_models,
            dataset_type=args.dataset_type,
            loss_function=args.loss_function,
            uncertainty_dropout_p=args.uncertainty_dropout_p,
            dropout_sampling_size=args.dropout_sampling_size,
            spectra_phase_mask=getattr(train_args, "spectra_phase_mask", None),
        )

    # Edge case if empty list of smiles is provided
    if len(test_data) == 0:
        preds = [None] * len(full_data)
        unc = [None] * len(full_data)
    else:
        preds, unc = predict_and_save(
            args=args,
            train_args=train_args,
            test_data=test_data,
            task_names=task_names,
            num_tasks=num_tasks,
            test_data_loader=test_data_loader,
            full_data=full_data,
            full_to_valid_indices=full_to_valid_indices,
            models=models,
            scalers=scalers,
            num_models=num_models,
            calibrator=calibrator,
            return_invalid_smiles=return_invalid_smiles,
            save_results = False,
        )

    if return_index_dict:
        preds_dict = {}
        unc_dict = {}
        for i in range(len(full_data)):
            if return_invalid_smiles:
                preds_dict[i] = preds[i]
                unc_dict[i] = unc[i]
            else:
                valid_index = full_to_valid_indices.get(i, None)
                if valid_index is not None:
                    preds_dict[i] = preds[valid_index]
                    unc_dict[i] = unc[valid_index]
        if return_uncertainty:
            return preds_dict, unc_dict
        else:
            return preds_dict
    else:
        if return_uncertainty:
            return preds, unc
        else:
            return preds
            

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
    cat_list,solv_list,reag_list = condition_key
    top3_indices = []
    for i in range(len(MPNN_out)):
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

    from Score import cal_acc
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path',type=str,default=args.test_path)
    parser.add_argument('--pred_path',type=str,default=args.save_path)
    args = parser.parse_args()
    print(n_list)
    cal_acc(args)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='reaction condition prediction')
    parser.add_argument('--test_path', type=str, default='./data/data_test.csv', help='path to test data')
    parser.add_argument('--model_path', type=str, default='./models', help='path to model')
    parser.add_argument('--label_path', type=str, default='./data/labels', help='path to condition keys')
    parser.add_argument('--library_path', type=str, default='./data/condition_library', help='path to classed conditions library')
    parser.add_argument('--save_path', type=str, default='./data/prediction', help='path to save condition prediction')
    args = parser.parse_args()
    Prediction(args)
