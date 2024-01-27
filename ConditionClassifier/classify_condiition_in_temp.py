import rdkit
from rdkit import Chem
import json
import csv
import gzip
import pandas as pd
from joblib import Parallel, delayed
from operator import add
from functools import reduce
from collections import Counter
import sys
import warnings
warnings.filterwarnings("ignore")

def Categorization_conditions(condition:str):
    category_dic = {'alkene':['[CX3:1]=[CX3:2]'],
                    'alkyne':['[CX2:1]#[CX2:2]'],
                    'alcohol':['[CX4:1][OX2H]'],
                    'ether':['[CX4:1][OX2][CX4:2]'],
                    'aldehyde':['[CX3H1:1][CX3H1](=[OX1])[#6:2]'],
                    'ketone':['[CX3H1:1][CX3H1](=[OX1])[#6:2]'],
                    'carboxylic acid':['[!O:1][CX3:2](=[OX1])[OX2H1:3]'],
                    'ester':['[CX3:1](=[OX1])[OX2H0:2]'],
                    'amide':['[CX3:1](=[OX1])[NX3H2:2]'],
                    'nitro':['[NX3:1](=[OX1])(=[OX1])[OX2:2]'],
                    'amine':['[NX3:1]','[n]'],
                    'halide':['[F;H0;X1]','[Cl;H0;X1]','[Br;H0;X1]', '[I;H0;X1]'],
                    'acid chloride':['[CX3:1](=[OX1])[Cl,Br,I:2]'],
                    'anhydride':['[CX3:1](=[OX1])[OX2:2][CX3:3](=[OX1])[OX2:4]'],
                    'nitrile':['[NX1-:1]#[CX2+0:2]'],
                    'aromatic':['a'],
                    'sulfone/sulfoxide':['[SX4](=[OX1])(=[OX1])','[SX3](=[OX1])'],
                    'phosphine' : ['[PX3]','[PX4]','[PX5]','[PX6]'],
                    'transition metal':['[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,#39,#40,#41,#42,#44,#45,#46,#47,#48,#72,#73,#74,#75,#76,#77,#78,#79,#80,#104,#105,#106,#107,#108,#109,#110,#111,#112]'],
                    'reducing metal':['[#3,#4,#11,#12,#13,#19,#20,#30,#26;+0]'],
                    'Main group metal':['[#13,#31,#49,#81,#50,#82,#83;+0,+1]'],
                    'Metal oxidizer':['[Cr+6]','[Mn+7]','[Mn+4]','[Ce+4]','[Pb+4]'],
                    'reductant':['[H-]','[BH4-]','[AlH4-]','[NaH]','[LiH]','[BH3]','[BH2]','[BH]','[AlH3]','[AlH2]','[AlH]'],
                    'acid':['[ClH1,BrH1,IH1:1]','O=S(=O)(O)O','O=C(O)C(F)(F)F','[CX3:1](=[OX1])[OX2H1:2]'],
                    'lewis acid':['[AlX3:1][F,Cl,Br,I,C:2]','[BX3:1][F,Cl,Br,I,H:2]','[Al+3]','[Ti+4]','[Zn+2]','[ZnX2:1][Cl,Br,I:2]','[Si+4]','[Fe+3]','[FeX3:1][Cl,Br,F:2]','[Ge+4]','[Sn+4]','[Ce+4]'],
                    'metal alkyl':['[CX4:1][Mg,Al,Zn,Li:2]'],
                    'silane':['[SiX4:1][#6:2]'],
                    'sulfide':['[CX4:1][SX2:2]'],
                    }
    out_list = []
    if condition == 'None':
        return []
    mol = Chem.MolFromSmiles(condition)
    for i in category_dic:
        for j in category_dic[i]:
            patt = Chem.MolFromSmarts(j)
            try:
                if mol.HasSubstructMatch(patt):
                    out_list.append(i)
                    continue
            except:
                pass
    if len(out_list) == 0:
        try:
            patt = Chem.MolFromSmarts('[+,-]')
            if mol.HasSubstructMatch(patt):
                out_list.append('ionic')
            else:
                patt = Chem.MolFromSmarts('[CX4:1]')
                if mol.HasSubstructMatch(patt):
                    out_list.append('alkane')
                else:
                    pass
        except:
            pass
    return list(set(out_list))

def get_condition_key(path:str):
    with open('%s/all_cat_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        cat_list_N = [row['cat'] for row in reader]

    with open('%s/all_solv_withN.csv'%path,'r') as f:
        reader = csv.DictReader(f)
        solv_list_N = [row['solv'] for row in reader]

    with open('%s/all_reag_withN.csv'%path,'r') as f:
            reader = csv.DictReader(f)
            reag_list_N = [row['reag'] for row in reader]
    
    return cat_list_N,solv_list_N,reag_list_N


def get_labeled_condition(test_list):
    out_list = []
    for i in test_list:
        out_list.append([i,[[Categorization_conditions(i[0])],
                            [Categorization_conditions(i[1]),Categorization_conditions(i[2])],
                            [Categorization_conditions(i[3]),Categorization_conditions(i[4]),Categorization_conditions(i[5])]],
        ])
    return out_list

def ismach(list1,list2,num:int = 2):
    solv_reag1 = list1[1]+list1[2]
    solv_reag2 = list2[1]+list2[2]
    if str(list2) != '[[], [], []]' and ((set(list1[0][0])&set(list2[0]) != set() or list1[0][0]==list2[0]) and (len(set(reduce(add,solv_reag1))&set(solv_reag2))>=num or (set(reduce(add,solv_reag1))&set(solv_reag2)==set(solv_reag2) and set(solv_reag2) !=set()))):
        return True
    return False

def how_mach(list1,list2):
    solv_reag1 = list1[1]+list1[2]
    solv_reag2 = list2[1]+list2[2]
    return len(set(list1[0][0])&set(list2[0]))/max(len(set(list2[0])),1)+len(set(reduce(add,solv_reag1))&set(solv_reag2))/max(len(set(solv_reag2)),1)

def update_class_label(condition_label, args):
    cat_labels, solv_labels, reag_labels = [], [], []
    for i in condition_label:
        cat_labels += i[0][0]
        solv_labels += reduce(add,i[1])
        reag_labels += reduce(add,i[2])
    cat_count = Counter(cat_labels)
    solv_count = Counter(solv_labels)
    reag_count = Counter(reag_labels)
    out = [[],[],[]]
    for i in cat_count:
        if cat_count[i] > args.Inclusion*len(condition_label):
            out[0].append(i)
    for i in solv_count:
        if solv_count[i] > args.Inclusion*len(condition_label):
            out[1].append(i)
    for i in reag_count:
        if reag_count[i] > args.Inclusion*len(condition_label):
            out[2].append(i)
    return out


def classify(args, condition_list, out_list: list = [],num:int = 2, update_label:bool = True):
    for i in range(len(condition_list)):
        if len(out_list) == 0:
            out_list.append({'class_label':None, 'conditions':[condition_list[i][0]], 'condition_label':[condition_list[i][1]]})
            out_list[0]['class_label'] = update_class_label([condition_list[i][1]], args)
        else:
            mach_list = []
            for j in range(len(out_list)):
                if ismach(condition_list[i][1],out_list[j]['class_label'],num):
                    mach_list.append([j,how_mach(condition_list[i][1],out_list[j]['class_label'])])
            if len(mach_list) == 0:
                out_list.append({'class_label':None, 'conditions':[condition_list[i][0]], 'condition_label':[condition_list[i][1]]})
                out_list[-1]['class_label'] = update_class_label([condition_list[i][1]], args)
            else:
                mach_list.sort(key=lambda x:x[1])
                out_list[mach_list[0][0]]['conditions'].append(condition_list[i][0])
                out_list[mach_list[0][0]]['condition_label'].append(condition_list[i][1])
                if update_label or len(out_list[mach_list[0][0]]['conditions'])<10:
                    out_list[mach_list[0][0]]['class_label'] = update_class_label(out_list[mach_list[0][0]]['condition_label'], args)

    return out_list

def merge_list(class_list):
    all_class = []
    for i in class_list:
        if len(all_class) == 0:
            all_class.append(i)
        else:
            for j in range(len(all_class)):
                set1 = {item for sublist in all_class[j]['class_label'] for item in sublist}
                set2 = {item for sublist in i['class_label'] for item in sublist}
                if set1 == set2:
                    if all_class[j]['class_label'][0] == set() and all_class[j]['class_label'][2] == set():
                        all_class[j]['class_label'] = i['class_label']
                    all_class[j]['conditions'] += i['conditions']
                    all_class[j]['condition_label'] += i['condition_label']
                    break
            else:
                all_class.append(i)
    return all_class


def get_unnecessary_condition(conditions):
    all_conditions = reduce(add,conditions)
    all_conditions = Counter(all_conditions)
    unnecessary_conditions = [i if (all_conditions[i] == 1 and all_conditions[i]<0.1*len(conditions)) else None for i in all_conditions]
    return [i for i in unnecessary_conditions if i != None]



def remove_unnecessary_condition(condition_list):
    for i in range(len(condition_list)):
        unnecessary_conditions = get_unnecessary_condition(condition_list[i]['conditions'])
        for j in range(len(condition_list[i]['conditions'])):
            condition = condition_list[i]['conditions'][j]
            for reagent in condition:
                if reagent in unnecessary_conditions:
                    condition[condition.index(reagent)] = 'None'
            new_label = get_labeled_condition([condition])[0][1]
            condition_list[i]['condition_label'][j] = new_label
        condition_list[i]['class_label'] = update_class_label(condition_list[i]['condition_label'])
    return condition_list

def reclassify(condition_list,args):
    classed_list = []
    reclass_list = []
    for i in range(len(condition_list)):
        if len(reduce(add,condition_list[i]['class_label']))>0:
            classed_list.append({'class_label':condition_list[i]['class_label'],'conditions':condition_list[i]['conditions'],'condition_label':condition_list[i]['condition_label']})
        else:
            reclass_list += condition_list[i]['conditions']
    reclass_list = get_labeled_condition(reclass_list)
    reclass_list = classify(args,reclass_list,classed_list,2,update_label = False)
    return reclass_list

def check(class_list):
    o1 = 0
    for i in class_list:
        o1 += len(i['conditions'])
    return o1

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
        try:
            cat = cat_list.index(condition[0])
        except:
            cat = 0
        try:
            solv1 = solv_list.index(condition[1])
        except:
            solv1 = 0
        try:
            solv2 = solv_list.index(condition[2])
        except:
            solv2 = 0
        try:
            reag1 = reag_list.index(condition[3])
        except:
            reag1 = 0
        try:
            reag2 = reag_list.index(condition[4])
        except:
            reag2 = 0
        try:
            reag3 = reag_list.index(condition[5])
        except:
            reag3 = 0
        out.append([cat,solv1,solv2,reag1,reag2,reag3])
    return out


def Classify_reaction_conditions(test_list,tem,smart,args):
    cat_list,solv_list,reag_list = get_condition_key('data')
    test_list = [eval(str(i)) for i in test_list]
    test_list = get_labeled_condition(test_list)
    test_list = classify(args,test_list,[],2)
    test_list = merge_list(test_list)
    test_list = reclassify(test_list,args)
    test_list = merge_list(test_list)
    out_list = []
    for i in range(len(test_list)):
        try:
            encoded_conditions = encode_condition(test_list[i]['conditions'],cat_list,solv_list,reag_list)
        except:
            print('___________________________________________________________')
        out_list.append({'tem':str(tem),'tem_smart':smart,'class_id':str(tem)+'_%s'%i,'class_label':test_list[i]['class_label'],'conditions':test_list[i]['conditions'],'encoded_conditions':encoded_conditions})
    return out_list



if __name__ == '__main__':
    '''
    condition_list = get_tem_condition('data')
    classed_condition_list = list(Parallel(n_jobs=-1)(delayed(Classify_reaction_conditions)(sorted(condition_list['conditions'][i]),condition_list['tem'][i],condition_list['tem_smart'][i]) for i in range(len(condition_list))))
    with open('data/classed_condition_list.json','w') as f:
        json.dump(classed_condition_list,f)
    '''
    pass


    
    