import csv
import pandas as pd
from rdkit import Chem
from sklearn.model_selection import train_test_split
import torch
import torch.utils.data as Data
from torch import nn
from torch.optim import SGD,Adam
from MLP_module import MLP_module0, MLP_module1, MLP_module2, MLP_module3,MLP_module5
csv.field_size_limit(500 * 1024 * 1024)
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
import numpy as np
import pandas as pd
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdChemReactions import RemoveMappingNumbersFromReactions
import time

with open('data/all_cat3.csv','r') as f:
    reader = csv.DictReader(f)
    for classes in reader:
        cat_list = list(classes.keys())[:866]

with open('data/all_solv3.csv','r') as f:
    reader = csv.DictReader(f)
    for classes in reader:
        solv_list = list(classes.keys())[:2702]

#None or Not None
cat_list.append('None')
solv_list.append('None')
cat_list.append('else')
solv_list.append('else')


with open('data/temp_condition.csv','r') as f:
    reader = csv.DictReader(f)
    dic_list = list(reader)

def get_one_hot(tem):
    blist = []
    for i in list(tem):
        alist = [0]*25065
        alist[int(i)] = 1
        blist.append(alist)
    tems = torch.tensor(blist, dtype=torch.float32)
    return tems

def get_top_1(outputs, tem, target, dlist = dic_list):
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


def get_top_3(outputs, tem, target, dlist = dic_list):
    out = []
    for i in range(outputs.size(0)):
        t_list = eval(dlist[int(tem[i].data)][target])
        Max_index = [-1]*3
        Max = [-1000,-1000,-1000]
        for j in range(outputs[i].size(0)):
            if outputs[i][j].data > Max[0] and j in t_list:
                Max[0] = outputs[i][j].data
                Max_index[0] = j
            elif outputs[i][j].data > Max[1] and j in t_list:
                Max[1] = outputs[i][j].data
                Max_index[1] = j
            elif outputs[i][j].data > Max[2] and j in t_list:
                Max[2] = outputs[i][j].data
                Max_index[2] = j
        out.append(Max_index)
    return torch.tensor(out)


def get_top_10(outputs, tem, target, dlist = dic_list):
    out = []
    for i in range(outputs.size(0)):
        t_list = eval(dlist[int(tem[i].data)][target])
        Max_index = [-1]*10 
        Max = [-1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000]
        for j in range(outputs[i].size(0)):
            if outputs[i][j].data > Max[0] and j in t_list:
                Max[0] = outputs[i][j].data
                Max_index[0] = j
            elif outputs[i][j].data > Max[1] and j in t_list:
                Max[1] = outputs[i][j].data
                Max_index[1] = j
            elif outputs[i][j].data > Max[2] and j in t_list:
                Max[2] = outputs[i][j].data
                Max_index[2] = j
            elif outputs[i][j].data > Max[3] and j in t_list:
                Max[3] = outputs[i][j].data
                Max_index[3] = j
            elif outputs[i][j].data > Max[4] and j in t_list:
                Max[4] = outputs[i][j].data
                Max_index[4] = j
            elif outputs[i][j].data > Max[5] and j in t_list:
                Max[5] = outputs[i][j].data
                Max_index[5] = j
            elif outputs[i][j].data > Max[6] and j in t_list:
                Max[6] = outputs[i][j].data
                Max_index[6] = j
            elif outputs[i][j].data > Max[7] and j in t_list:
                Max[7] = outputs[i][j].data
                Max_index[7] = j
            elif outputs[i][j].data > Max[8] and j in t_list:
                Max[8] = outputs[i][j].data
                Max_index[8] = j
            elif outputs[i][j].data > Max[9] and j in t_list:
                Max[9] = outputs[i][j].data
                Max_index[9] = j
        out.append(Max_index)
    return torch.tensor(out)


def train_module3(model,target, train_loader,test_loader,loss_function = nn.CrossEntropyLoss().to(device=torch.device("cuda" if torch.cuda.is_available() else "cpu")),Ir = 0.0001,epochs = 5):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    optimizer = torch.optim.Adam(model.parameters(),lr=Ir)
    for epoch in range(epochs):
        start = time.time()
        running_loss = 0.0
        for step,(b_x,b_p,b_tem,b_tar) in enumerate(train_loader):
            b_x = b_x.to(device)
            b_p = b_p.to(device)
            b_tem = b_tem.to(device)
            b_tar = b_tar.to(device)
            optimizer.zero_grad()
            b_t=  get_one_hot(b_tem)
            b_t = b_t.to(device)
            out = model(b_x,b_p,b_t)
            loss = loss_function(out,b_tar)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
            if step % 900 == 899:
                print('[%d, %5d] loss: %.3f' % (epoch + 1, step + 1, running_loss / 900))
                running_loss = 0.0
        end = time.time()
        test_module3(model,target,test_loader,use_all=False)
        print('time:',end - start)
        torch.save(model,'%s_model3_with_None.pt'%target)

def test_module3(model,target,test_loader,use_all = True):
    correct = 0
    total = 0
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    with torch.no_grad(): 
        for step,(b_x,b_p,b_tem,b_tar) in enumerate(test_loader):
            b_x = b_x.to(device)
            b_p = b_p.to(device)
            b_tem = b_tem.to(device)
            b_tar = b_tar.to(device)
            b_t=  get_one_hot(b_tem)
            b_t = b_t.to(device)
            outputs = model(b_x,b_p,b_t)
            _,predicted = torch.max(outputs.data,dim = 1) 
            total += b_tar.size(0)
            correct += (predicted == b_tar).sum().item()
            if use_all == False and step >=100:
                break
    print('Accuracy on test set: %.4f' % (correct / total))
    return correct / total

def topk_acc3(model,target,test_loader):
    #device = torch.device("")
    correct_0 = 0
    correct_1 = 0
    correct_3 = 0
    correct_10 = 0
    total = 0
    with torch.no_grad(): 
        for step,(b_x,b_p,b_tem,b_tar) in enumerate(test_loader):
            '''
            b_x = b_x.to(device)
            b_p = b_p.to(device)
            b_tem = b_tem.to(device)
            b_tar = b_tar.to(device)
            '''
            b_t=  get_one_hot(b_tem)
            #b_t = b_t.to(device)
            outputs = model(b_x,b_p,b_t)
            #_, top_10 = torch.topk(outputs.data,10) 
            #_, top_3 = torch.topk(outputs.data,3)
            _, pre_1 = torch.max(outputs.data,dim = 1)
            top_1 = get_top_1(outputs.data,b_tem,target)
            top_3 = get_top_3(outputs.data,b_tem,target)
            top_10 = get_top_10(outputs.data,b_tem,target)
            total += b_tar.size(0)
            correct_0 += (pre_1 == b_tar).sum().item()
            correct_1 += (top_1 == b_tar).sum().item()
            for i in range(b_tar.size(0)):
                if b_tar[i] in top_3[i]:
                    correct_10 += 1
                    correct_3 += 1
                elif b_tar[i] in top_10[i]:
                    correct_10 += 1
            
    print('MLP3-%s'%target)
    print('pre-1 acc:',correct_0/total)
    print('top-1 acc: ',correct_1/total)
    print('top-3 acc: ',correct_3/total)
    print('top-10 acc: ',correct_10/total)


def test_MLP_module_withoutNone(target):
    data = pd.read_json('data/MLP_%s_data.json.gz'%target,compression='gzip')
    X_train, X_test, y_train, y_test = train_test_split(data[['rfpgen','pfpgen','tem']], data[[target]], test_size=0.1,random_state= 2)
    X_train_tensor0 = torch.tensor(list(X_train['rfpgen']), dtype=torch.float32)
    X_train_tensor1 = torch.tensor(list(X_train['pfpgen']), dtype=torch.float32)
    X_train_tensor2 = torch.tensor(list(X_train['tem']), dtype=torch.float32)
    train_tensor = torch.tensor(list(y_train[target]), dtype=torch.int64)
    X_test_tensor0 = torch.tensor(list(X_test['rfpgen']), dtype=torch.float32)
    X_test_tensor1 = torch.tensor(list(X_test['pfpgen']), dtype=torch.float32)
    X_test_tensor2 = torch.tensor(list(X_test['tem']), dtype=torch.float32)
    test_tensor = torch.tensor(list(y_test[target]), dtype=torch.int64)
    train_dataset = Data.TensorDataset(X_train_tensor0, X_train_tensor1,X_train_tensor2,train_tensor)  
    test_dataset = Data.TensorDataset(X_test_tensor0,X_test_tensor1, X_test_tensor2,test_tensor)
    train_loader = Data.DataLoader(dataset=train_dataset,batch_size=100,shuffle=True,num_workers=0)
    test_loader = Data.DataLoader(dataset=test_dataset,batch_size=100,shuffle=True,num_workers=0)
    print('data get!')
    '''
    if target == 'cat':
        model = MLP_module3(866)
        epochs = 40
    if target == 'solv':
        model = MLP_module3(2702)
        epochs = 5
    '''
    #device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = torch.load('%s_model3_without_None.pt'%target)
    epochs = 5
    #train_module3(model,target, train_loader,test_loader,loss_function = nn.CrossEntropyLoss(),Ir = 0.0001,epochs =epochs)
    print('------------------------------------')
    topk_acc3(model,target,test_loader)
    print('------------------------------------')

#------------------------------------------------------------------

def test_MLP_module_withNone(target):
    
    data = pd.read_json('data/MLP_all_data.json.gz',compression='gzip')
    X_train, X_test, y_train, y_test = train_test_split(data[['rfpgen','pfpgen','tem']], data[[target]], test_size=0.1,random_state= 2)
    X_train_tensor0 = torch.tensor(list(X_train['rfpgen']), dtype=torch.float32)
    X_train_tensor1 = torch.tensor(list(X_train['pfpgen']), dtype=torch.float32)
    X_train_tensor2 = torch.tensor(list(X_train['tem']), dtype=torch.float32)
    train_tensor = torch.tensor(list(y_train[target]), dtype=torch.int64)
    X_test_tensor0 = torch.tensor(list(X_test['rfpgen']), dtype=torch.float32)
    X_test_tensor1 = torch.tensor(list(X_test['pfpgen']), dtype=torch.float32)
    X_test_tensor2 = torch.tensor(list(X_test['tem']), dtype=torch.float32)
    test_tensor = torch.tensor(list(y_test[target]), dtype=torch.int64)
    train_dataset = Data.TensorDataset(X_train_tensor0, X_train_tensor1,X_train_tensor2,train_tensor)  
    test_dataset = Data.TensorDataset(X_test_tensor0,X_test_tensor1, X_test_tensor2,test_tensor)
    train_loader = Data.DataLoader(dataset=train_dataset,batch_size=100,shuffle=True,num_workers=0)
    test_loader = Data.DataLoader(dataset=test_dataset,batch_size=100,shuffle=True,num_workers=0)
    print('train data get!')
    
    if target == 'cat':
        model = MLP_module3(868)
        epochs = 4
    if target == 'solv':
        model = MLP_module3(2704)
        epochs = 4
    
    
    
    model = torch.load('%s_model3_with_None.pt'%target)
    #train_module3(model,target, train_loader,test_loader,loss_function = nn.CrossEntropyLoss(),Ir = 0.0001,epochs =epochs)


    data = pd.read_json('data/MLP_%s_data.json.gz'%target,compression='gzip')
    X_train, X_test, y_train, y_test = train_test_split(data[['rfpgen','pfpgen','tem']], data[[target]], test_size=0.1,random_state= 2)
    X_train_tensor0 = torch.tensor(list(X_train['rfpgen']), dtype=torch.float32)
    X_train_tensor1 = torch.tensor(list(X_train['pfpgen']), dtype=torch.float32)
    X_train_tensor2 = torch.tensor(list(X_train['tem']), dtype=torch.float32)
    train_tensor = torch.tensor(list(y_train[target]), dtype=torch.int64)
    X_test_tensor0 = torch.tensor(list(X_test['rfpgen']), dtype=torch.float32)
    X_test_tensor1 = torch.tensor(list(X_test['pfpgen']), dtype=torch.float32)
    X_test_tensor2 = torch.tensor(list(X_test['tem']), dtype=torch.float32)
    test_tensor = torch.tensor(list(y_test[target]), dtype=torch.int64)
    train_dataset = Data.TensorDataset(X_train_tensor0, X_train_tensor1,X_train_tensor2,train_tensor)  
    test_dataset = Data.TensorDataset(X_test_tensor0,X_test_tensor1, X_test_tensor2,test_tensor)
    train_loader = Data.DataLoader(dataset=train_dataset,batch_size=100,shuffle=True,num_workers=0)
    test_loader = Data.DataLoader(dataset=test_dataset,batch_size=100,shuffle=True,num_workers=0)
    print('test data get!')

    print('------------------------------------')
    topk_acc3(model,target,test_loader)
    print('------------------------------------')

#------------------------------------------------------------------

def popular_0(target):
    if target == 'cat':
        return [0,1,2,3,4,5,6,7,8,9]
    if target == 'solv':
        return [0,1,2,3,4,5,6,7,8,9]

def popular_1(x,target):
    pdata = pd.read_csv('./data/temp_popular_condition.csv')
    if target == 'cat':
        return eval(pdata['cat'][x])
    if target == 'solv':
        return eval(pdata['solv'][x])

    
def try_popular_model(target):
    data = pd.read_json('data/MLP_%s_data.json.gz'%target,compression='gzip')
    X_train, X_test, y_train, y_test = train_test_split(data[['rfpgen','pfpgen','tem']], data[[target]], test_size=0.1,random_state= 2)
    p0_top1 = 0
    p0_top3 = 0
    p0_top10 = 0
    p1_top1 = 0
    p1_top3 = 0
    p1_top10 = 0
    total = len(X_test['tem'])
    ys = list(y_test[target])
    xs = list(X_test['tem'])
    for i in range(len(X_test['tem'])):
        if i % (len(X_test['tem'])//100) == 0:
            print((i*100)/len(X_test['tem']),'%')
        y = int(ys[i])
        out0 = popular_0(target)
        out1 = popular_1(xs[i],target)
        if y == out0[0]:
            p0_top1 += 1
            p0_top3 += 1
            p0_top10 += 1
        elif y in out0[:3]:
            p0_top3 += 1
            p0_top10 += 1
        elif y in out0:
            p0_top10 += 1
        if y == out1[0]:
            p1_top1 += 1
            p1_top3 += 1
            p1_top10 += 1
        elif y in out1[:3]:
            p1_top3 += 1
            p1_top10 += 1
        elif y in out1:
            p1_top10 += 1
    
    print('-----------------------------------------------')
    print('popular-0-%s'%target)
    print('top-1 acc: ',p0_top1/total)
    print('top-3 acc: ',p0_top3/total)
    print('top-10 acc: ',p0_top10/total)
    print('-----------------------------------------------')
    print('popular-1-%s'%target)
    print('top-1 acc: ',p1_top1/total)
    print('top-3 acc: ',p1_top3/total)
    print('top-10 acc: ',p1_top10/total)
    print('-----------------------------------------------')



test_MLP_module_withoutNone('cat')
test_MLP_module_withoutNone('solv')

#try_popular_model('cat')
#try_popular_model('solv')


