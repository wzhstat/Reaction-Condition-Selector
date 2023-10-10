import csv
import pandas as pd
from rdkit import Chem
from sklearn.model_selection import train_test_split
import torch
import torch.utils.data as Data
from torch import nn
from torch.optim import SGD,Adam
from .MLPModel import nnModel0, nnModel1, nnModel2
from joblib import Parallel, delayed
csv.field_size_limit(500 * 1024 * 1024)
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
import numpy as np
import pandas as pd
import time
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdChemReactions import RemoveMappingNumbersFromReactions
from collections import defaultdict

def get_one_hot_tem(tem, teml):
    """
    Get the one-hot encoding of the template.

    Args:
        tem (torch.Tensor): The template.
        teml (int): The length of the template.

    Returns:
        torch.Tensor: The one-hot encoding of the template.
    """
    # Create a list of zeros with length teml for each input in the template.
    blist = [[0] * teml for _ in range(len(tem))]

    # Set the position of the input value to 1 in each list.
    for i, t in enumerate(tem):
        blist[i][int(t)] = 1

    # Convert the list to a tensor.
    return torch.tensor(blist, dtype=torch.float32)


def get_top_1(add_filter, tem, outputs, target, dlist):
    """
    Get the top 1 prediction for each input.

    Args:
        add_filter (bool): Whether to add a filter.
        tem (torch.Tensor): The template.
        outputs (torch.Tensor): The outputs.
        target (torch.Tensor): The target.
        dlist (List[str]): The list of filters.

    Returns:
        torch.Tensor: The top 1 prediction.
    """
    if add_filter:
        out = []
        for i in range(outputs.size(0)):
            t_list = eval(dlist[int(tem[i].data)][target])
            Max_index = None
            Max = -1000
            for j in range(outputs[i].size(0)):
                if outputs[i][j].data > Max and j in t_list:  
                    Max = outputs[i][j].data
                    Max_index = j
            if Max_index == None:
                Max_index = t_list[0]
            out.append(Max_index)
        return torch.tensor(out)
    else:
        _, predicted = torch.max(outputs.data, dim = 1) 
        return predicted


def get_rxnfp(reaction):
    '''
    Get the fingerprint of the reaction.

    Args:
        reaction (tuple):reactant and product smiles
    
    Returns:
        tuple: reactant fingerprint, product fingerprint, reaction fingerprint
    '''
    try:
        # Get the reactants and products from the reaction.
        (reactant,product) = reaction
        rm = Chem.MolFromSmiles(reactant)
        pm = Chem.MolFromSmiles(product)

        # Get the Morgan fingerprint for each molecule.
        info = {}
        rfpgen= np.array(AllChem.GetMorganFingerprintAsBitVect(rm, useChirality=True, radius=2, nBits = 512, bitInfo=info))
        pfpgen= np.array(AllChem.GetMorganFingerprintAsBitVect(pm, useChirality=True, radius=2, nBits = 512, bitInfo=info))
        
        # Calculate the reaction fingerprint by subtracting the product fingerprint from the reactant fingerprint.
        rxnfp = pfpgen-rfpgen
        return (rfpgen,pfpgen,rxnfp)
    except Exception as e:
        print(e)
        return (None,None,None)



def get_conditionfp(condition):
    """
    Get the condition fingerprint for a given condition.

    Args:
        condition (str): The condition smiles

    Returns:
        numpy.ndarray: The condition fingerprint.
    """
    try:
        condition = Chem.MolFromSmiles(condition)
        info = {}
        confp = np.array(AllChem.GetMorganFingerprintAsBitVect(condition, useChirality=True, radius=2, nBits = 512, bitInfo=info))
        return confp
    except Exception as e:
        print(e)
        return np.array([0]*512)

def sum_tem_condition(withN, path, file_name):
    """
    Summarize the condition of each template in the given file.

    Args:
        withN (bool): Whether to include nitrogen in the condition.
        path (str): The path to the file.
        file_name (str): The name of the file.

    Returns:
        dict: A dictionary containing the condition of each template.
    """
    # Set the suffix for the file name based on whether nitrogen is included in the condition.
    if withN:
        N = 'withN'
    else:
        N = 'withoutN'

    # Define the list of conditions to count.
    count_conditions = ['cat', 'solv', 'reag']

    # Loop through each condition and read the corresponding file.
    for count_condition in count_conditions:
        with open('%s/all_%s_%s.csv' % (path, count_condition, N), 'r') as f:
            reader = csv.DictReader(f)
            for classes in reader:
                # Get the list of classes for each condition.
                if count_condition == 'cat':
                    cat_list = list(classes.keys())
                elif count_condition == 'solv':
                    solv_list = list(classes.keys())
                else:
                    reag_list = list(classes.keys())

    # Read the input file.
    data = pd.read_csv('%s/%s.csv' % (path, file_name))
    l = len(data['template'])
    all_data = {}

    # Loop through each template in the input file.
    for i in range(l):
        # Print the progress every 1%.
        if i % (l // 100) == 0:
            print('%d%%' % (i / (l // 100)))

        # Get the template.
        tem = data['template'][i]

        # If the template is not in the dictionary, add an empty set each condition.
        conditions = ['cat','solv','reag0','reag1','reag2']
        if tem not in all_data:
            all_data[tem] = {c: set() for c in ['cat','solv','reag0','reag1','reag2']}
        # Add the condition to the set for each condition.
        for condition in conditions:
            if condition in ['reag0','reag1','reag2']:
                if data[condition][i] in reag_list:
                    all_data[tem][condition].add(data[condition][i])
            else:
                if data[condition][i] in eval(condition+'_list'):
                    all_data[tem][condition].add(data[condition][i])
        all_data = pd.DataFrame(all_data).T
        all_data.to_csv('%s/temp_condition.csv' % path)
    return all_data


def generate_train_data(
        inputs: str,
        path: str,
        file_name: str,
        withN: bool,
        target: str
        ):
    '''
    Generate training data from the specified target.
    Args:
        inputs (str): The input to the model.
        path (str): The path to the data file.
        file_name (str): The name of the data file.
        withN (bool): Whether to include None in the condition.
        target (str): The target to predict.
    Returns:
        tuple: The training data, the number of classes, and the number of templates.
    '''

    # Load the data file.
    data = pd.read_csv('%s/%s.csv'%(path,file_name))

    # Set the suffix for the file name based on whether None is included in the condition.
    if withN:
        file_name1 = "withN"
    else:
        file_name1 = "withoutN"

    if target in ['cat','solv']:
        with open('%s/all_%s_%s.csv'%(path,target,file_name1),'r') as f:
            reader = csv.DictReader(f)
            for classes in reader:
                target_list = list(classes.keys())
    elif target in ['reag0','reag1','reag2','reag3']:
        with open('%s/all_reag_%s.csv'%(path,file_name1),'r') as f:
            reader = csv.DictReader(f)
            for classes in reader:
                target_list = list(classes.keys())
    else:
        raise KeyError("target must be 'cat','solv','reag0','reag1','reag2','reag3'")

    print("n %s:"%target,len(target_list))

    # Generate the reaction fingerprints.
    rxnfps = Parallel(n_jobs=-1, verbose=4)(delayed(get_rxnfp)(reaction) for reaction in list(data[['reactants','products']].apply(tuple, axis=1)))
    t_data = []
    max_tem = 0

    # Generate the reaction fingerprints.
    conditions = ['cat','solv','reag0','reag1','reag2']
    conditionfps = {"catfp":[],"solvfp":[],"reag0fp":[],"reag1fp":[],"reag2fp":[]}
    for condition in conditions:
        if condition in inputs:
            conditionfps[condition+'fp'] = Parallel(n_jobs=-1, verbose=4)(delayed(get_conditionfp)(reaction) for reaction in data[condition])

    for i in range(len(rxnfps)):
        dic = {}
        if data[target][i] not in target_list:
            continue
        if rxnfps[i][0] is None:
            continue
        dic['input'] = np.concatenate((rxnfps[i][0],rxnfps[i][1]))
        if 'rxnfp' in inputs:
            dic['input'] =np.concatenate((dic['input'],rxnfps[i][2]))
        for condition in conditions:
            if condition in inputs:
                dic['input'] = np.concatenate((dic['input'],conditionfps[condition+'fp'][i]))
        dic['tem'] = data['template'][i]
        if data['template'][i] > max_tem:
            max_tem = data['template'][i]
        dic[target] = target_list.index(data[target][i])
        t_data.append(dic)
    print('n template:',max_tem+1)
    t_data = pd.DataFrame(t_data)
    return t_data,len(target_list),max_tem+1




def train_model(
        model: nn.Module,
        target: str, 
        train_loader: Data.DataLoader,
        test_loader:  Data.DataLoader,
        add_filter: bool,
        loss_function,
        Ir: float,
        epochs: int,
        withN: bool,
        dic_list:list):
    '''
    Train the reaction classification model.
    
    Args:
        model (nn.Module): The model to train.
        target (str): The target to predict.
        train_loader (Data.DataLoader): The training data.
        test_loader (Data.DataLoader): The test data.
        add_filter (bool): Whether to add a filter.
        loss_function (nn.Module): The loss function.
        Ir (int): The learning rate.
        epochs (int): The number of epochs.
        withN (bool): Whether to include None in the condition.
        dic_list (list): The list of filters.

    '''
    optimizer = torch.optim.Adam(model.parameters(),lr=Ir)
    for epoch in range(epochs):
        # Set the model to training mode.
        model.train()
        running_loss = 0.0

        # Loop through the training data in batches.
        for step,data in enumerate(train_loader):
            # Zero the gradients.
            optimizer.zero_grad()

            # Get inputs
            b_t = data[-1]

            # Forward pass.
            out = model(data[0])
            loss = loss_function(out,b_t)

            # Backward pass and optimization.
            loss.backward()
            optimizer.step()
            running_loss += loss.item()

            # Print the and validation accuracy.
            if step % 900 == 899:
                print('[%d, %5d] loss: %.3f' % (epoch + 1, step + 1, running_loss / 900))
                running_loss = 0.0
        acc = test_model(model,target,test_loader,add_filter,dic_list,use_all=False)

    # Save the trained model.
    if withN:
        torch.save(model,'models/%s_model_withN.pt'%target)
    else:
        torch.save(model,'models/%s_model_withoutN.pt'%target)
    return acc

def test_model(
        model: nn.Module,
        target: str, 
        test_loader: Data.DataLoader,
        add_filter: bool,
        dic_list: list,
        use_all: bool= True):
    '''
    Test the reaction classification model.

    Args:
        model (nn.Module): The trained model.
        target (str): The target reaction type.
        test_loader (DataLoader): The test data loader.
        add_filter (bool): Whether to use additional filters.
        dic_list (list): The list of dictionaries containing the filter information.
        use_all (bool): Whether to use all the test data.
    '''
    # Initialize the counters.
    correct = 0
    total = 0
    # Set the model to evaluation mode.
    with torch.no_grad(): 
        # Loop through the test data.
        for step,data in enumerate(test_loader):

            # Get inputs
            b_t = data[-1]

            # Get the predicted labels.
            outputs = model(data[0])
            predicted = get_top_1(add_filter,data[1],outputs.data,target,dic_list) 

            # Update the counters.
            total += b_t.size(0)
            correct += (predicted == b_t).sum().item()
            if use_all == False and step >=100:
                break
    # Print the and validation accuracy.
    print('Accuracy on test set: %.4f' % (correct / total))
    return correct / total


def topk_acc(
    model: nn.Module,
    test_loader: Data.DataLoader,
    k: int):

    '''
    This function is used to calculate the topk accuracy.

    Args:
        model (nn.Module): The trained model.
        test_loader (DataLoader): The test data loader.
        k (int): The value of k.

    '''
    # Initialize the counters.
    correct = 0
    total = 0
    # Set the model to evaluation mode.
    with torch.no_grad():

        # Loop through the test data.
        for step,data in enumerate(test_loader):

            # Get inputs
            b_t = data[-1]

            # Get the predicted labels.
            outputs = model(data[0])
            _, predicted = torch.topk(outputs.data, k = k, dim = 1)

            # Update the counters.
            total += b_t.size(0)
            for i in range(len(predicted)):
                if b_t[i] in predicted[i]:
                    correct += 1

    # Print the and validation accuracy.
    print('Top%d acc: %.4f' % (k,correct / total))
    return correct / total


def train_model_withT(
        model: nn.Module,
        target: str,
        train_loader: Data.DataLoader,
        test_loader:  Data.DataLoader,
        teml: int,
        add_filter: bool,
        loss_function,
        Ir:float,
        epochs: int,
        withN: bool,
        dic_list:list
        ):
    '''
    This function is used to train the model with template.

    Args:
        model (nn.Module): The model to train.
        target (str): The target to predict.
        train_loader (Data.DataLoader): The training data.
        test_loader (Data.DataLoader): The test data.
        teml (int): The length of the template.
        add_filter (bool): Whether to add a filter.
        loss_function (nn.Module): The loss function.
        Ir (int): The learning rate.
        epochs (int): The number of epochs.
        withN (bool): Whether to include None in the condition.
        dic_list (list): The list of filters.
    '''
    optimizer = torch.optim.Adam(model.parameters(),lr=Ir)
    for epoch in range(epochs):
        running_loss = 0.0
        for step,data in enumerate(train_loader):
            optimizer.zero_grad()
            b_t = data[-1]
            b_tem = get_one_hot_tem(data[1],teml)
            out = model((data[0],b_tem))
            loss = loss_function(out,b_t)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
            if step % 900 == 899:
                print('[%d, %5d] loss: %.3f' % (epoch + 1, step + 1, running_loss / 900))
                running_loss = 0.0
        acc = test_model_withT(model,target,test_loader,teml,add_filter,dic_list,use_all=False)
    if withN:
        torch.save(model,'models/%s_model_withN.pt'%target)
    else:
        torch.save(model,'models/%s_model_withoutN.pt'%target)
    return acc



def test_model_withT(model,target, test_loader,teml,add_filter,dic_list,use_all = True):
    '''
    This function is used to test the model with template.

    Args:
        model (nn.Module): The trained model.
        target (str): The target reaction type.
        test_loader (DataLoader): The test data loader.
        teml (int): The length of the template.
        add_filter (bool): Whether to use additional filters.
        dic_list (list): The list of dictionaries containing the filter information.
        use_all (bool): Whether to use all the test data.
    '''
    correct = 0
    total = 0
    with torch.no_grad(): 
        for step,data in enumerate(test_loader):
            b_t = data[-1]
            b_tem = get_one_hot_tem(data[1],teml)
            outputs = model((data[0],b_tem))
            predicted = get_top_1(add_filter,data[1],outputs.data,target,dic_list) 
            total += b_t.size(0)
            correct += (predicted == b_t).sum().item()
            if use_all == False and step >=100:
                break
    print('Accuracy on test set: %.4f' % (correct / total))
    return correct / total


def topk_acc_withT(teml,model,test_loader,k):
    '''
    This function is used to calculate the topk accuracy with template.

    Args:
        model (nn.Module): The trained model.
        test_loader (DataLoader): The test data loader.
        k (int): The value of k.
    '''
    correct = 0
    total = 0
    with torch.no_grad():
        for step, data in enumerate(test_loader):
            b_t = data[-1]
            b_tem = get_one_hot_tem(data[1],teml)
            outputs = model((data[0],b_tem))
            _, predicted = torch.topk(outputs.data, k = k, dim = 1)
            total += b_t.size(0)
            for i in range(len(predicted)):
                if b_t[i] in predicted[i]:
                    correct += 1
    print('Top%d acc: %.4f' % (k,correct / total))
    return correct / total


def train(
        inputs: str,
        Model: nn.Module,
        path: str, 
        file_name: str, 
        withN: bool, 
        target: str, 
        epochs: int, 
        n1: int, 
        n2: int, 
        Ir: float, 
        batch_size: int, 
        add_filter: bool = False, 
        loss_function = nn.CrossEntropyLoss()):
    
    '''
    Train the reaction classification model.

    Args:
        inputs (str): The path to the input data.
        Model (nn.Module): The PyTorch model to train.
        path (str): The path to save the trained model.
        file_name (str): The name of the file to save the trained model.
        withN (bool): Whether to include nitrogen in the condition.
        target (str): The target reaction type.
        epochs (int): The number of epochs to train for.
        n1 (int): The size of the first hidden layer.
        n2 (int): The size of the second hidden layer.
        Ir (float): The initial learning rate.
        batch_size (int): The batch size.
        add_filter (bool): Whether to use additional filters.
        loss_function (nn.Module): The loss function to use.

    '''
    print('start to get train data')

    # Get the training data.
    data,targetl,teml = generate_train_data(inputs,path, file_name, withN, target)

    # Split the data into training and test sets.
    input1 = inputs.split('+')

    # Get the training and test data.
    if add_filter:
        with open('data/temp_condition.csv','r') as f:
            reader = csv.DictReader(f)
            dic_list = list(reader)
    else:
        dic_list = None
    X_train, X_test, y_train, y_test = train_test_split(data[['input','tem']], data[target], test_size=0.1)
    X_train_tensor0 = torch.tensor(list(X_train['input']), dtype=torch.float32)
    X_train_tensor2 = torch.tensor(list(X_train['tem']), dtype=torch.float32)
    y_train_tensor = torch.tensor(list(y_train), dtype=torch.int64)
    X_test_tensor0 = torch.tensor(list(X_test['input']), dtype=torch.float32)
    X_test_tensor2 = torch.tensor(list(X_test['tem']), dtype=torch.float32)
    y_test_tensor = torch.tensor(list(y_test), dtype=torch.int64)
    train_dataset = Data.TensorDataset(X_train_tensor0, X_train_tensor2, y_train_tensor)
    train_loader = Data.DataLoader(dataset=train_dataset, batch_size=batch_size, shuffle=True)
    test_dataset = Data.TensorDataset(X_test_tensor0, X_test_tensor2, y_test_tensor)
    test_loader = Data.DataLoader(dataset=test_dataset, batch_size=batch_size, shuffle=True)
    n0 = len(input1)*512
    print('n0:',n0)
    print('get data done')

    # Train the model.
    if Model in [nnModel0,nnModel1]:
        model = Model(targetl,n0,n1,n2)
        acc= train_model(model,target, train_loader,test_loader,add_filter = add_filter,loss_function = loss_function,Ir = Ir,epochs = epochs,withN = withN,dic_list=dic_list)
        print('------------------------------------')
        acc3 = topk_acc(model,test_loader,k=3)
        print('------------------------------------')
        acc10 = topk_acc(model,test_loader,k=10)
        print('------------------------------------')
        outdic = {'acc':acc,'acc3':acc3,'acc10':acc10}
        outdic = pd.DataFrame(outdic,index=[0])
        outdic.to_csv('models/%s_%s_out.csv'%(target,file_name))


    elif Model == nnModel2:
        model = Model(targetl,n0,n1,n2)
        acc= train_model_withT(model,target, train_loader,test_loader,teml,add_filter = add_filter,loss_function = loss_function,Ir = Ir,epochs = epochs,withN = withN,dic_list=dic_list)
        print('------------------------------------')
        acc3 = topk_acc_withT(model,test_loader,k=3)
        print('------------------------------------')
        acc10 = topk_acc_withT(model,test_loader,k=10)
        print('------------------------------------')
        outdic = {'acc':acc,'acc3':acc3,'acc10':acc10}
        outdic = pd.DataFrame(outdic,index=[0])
        outdic.to_csv('models/%s_%s_out.csv'%(target,file_name))
    
