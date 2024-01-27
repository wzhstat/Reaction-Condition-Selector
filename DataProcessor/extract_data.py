import os
from xml.dom import minidom
from xml.dom.minidom import parse
import csv
import gzip
import json
import pandas as pd
from rdkit import Chem
import sys
from rdchiral import template_extractor
from rdkit.Chem.rdChemReactions import RemoveMappingNumbersFromReactions
csv.field_size_limit(500 * 1024 * 1024)

'''
def extract(reaction):
    try:
        return template_extractor.extract_from_reaction(reaction)
    except KeyboardInterrupt:
        print('Interrupted')
        raise KeyboardInterrupt
    except Exception as e:
        print(e)
        return {'reaction_id': reaction['_id']}
'''

def open_csv(path:str):
    '''
    open csv file and get all cat,solv,reag
    args:
        path: csv file path
    '''
    
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


def remove_mapping(a):
    """
    Remove the mapping from Smarts.

    Args:
        a (str): A mol in Smarts.
    """
    a = Chem.MolFromSmiles(a)
    for atom in a.GetAtoms():
        atom.SetAtomMapNum(0)
    a = Chem.MolToSmiles(a)
    return a

def Cleam_merge(alist):
    '''
    Clean and merge a list of dictionaries.

    Args:
        alist (list): A list of dictionaries.

    Returns:
        A merged list.
    '''

    #Get the list of reactants and reagents.
    newlist = eval(str(alist))
    outlist = []
    charge_list = []

    #Clean and merge the list.
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
            if "+" in s or "-" in s:
                charge_list.append(s)
            else:
                outlist.append(s)
        charge_list.sort()
        outlist.sort()
    return charge_list+outlist


def check_charge(rlist,rglist):
    '''
    This function is used to check the charge of the reactants and reagents.

    Args:
        rlist (list): A list of reactants.
        rglist (list): A list of reagents.
    '''

    if len(rglist) == 0:
        return rlist,'None'
    else:
        charge = 0
        for i in rglist:
            for j in range(len(i)):
                if i[j] == '+':
                    if len(i) > j and i[j+1].isdigit():
                        charge += int(i[j+1])
                    else:
                        charge += 1
                elif i[j] == '-':
                    if len(i) > j and i[j+1].isdigit():
                        charge -= int(i[j+1])
                    else:
                        charge -= 1
        if charge == 0:
            return rlist,rglist
        else:
            nrglist = rglist.copy()
            if charge > 0:
                for i in rglist:
                    icharge = 0
                    for j in range(len(i)):
                        if i[j] == '+':
                            if len(i) > j and i[j+1].isdigit():
                                icharge += int(i[j+1])
                            else:
                                icharge += 1
                    if icharge > 0:
                        nrglist.remove(i)
                        rlist.append(i)
                        charge -= icharge
                    if charge == 0:
                        break
            else:
                for i in rglist:
                    icharge = 0
                    for j in range(len(i)):
                        if i[j] == '-':
                            if len(i) > j and i[j+1].isdigit():
                                icharge -= int(i[j+1])
                            else:
                                icharge -= 1
                    if icharge < 0:
                        nrglist.remove(i)
                        rlist.append(i)
                        charge += icharge
                    if charge == 0:
                        break
        return rlist,nrglist



def get_datas(args):
    out_path = "%s/%s.csv"%(args.save_path,args.data_name)
    load_path = args.data_path
    '''
    This function is used to get the data from the xml files and save them into a csv file.
    load_path: the path of the data
    save_path: the path of the csv file

    Args:
        data_name (str): The name of the dataset.
        load_path (str): The path of the data.
        save_path (str): The path to save the results.

    '''
    header = ['_id','source','reaction','products','reactants','reagent','catalyst','solvent']
    n = 1
    with open("%s/%s_all.csv"%(args.save_path,args.data_name), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)

        for filename in os.listdir(load_path):
            alllist = []
            for filename1 in os.listdir("%s/%s"%(load_path,filename)):
                doc = parse("%s/%s/%s"%(load_path,filename,filename1))
                reactionList = doc.documentElement
                reactions = reactionList.getElementsByTagName("reaction")
                for reaction in reactions:
                    alist = [] # list of a reaction
                    plist = [] # list of products
                    rlist = [] # list of reactants
                    rglist = [] # list of reactants
                    source = reaction.getElementsByTagName('dl:source')[0].getElementsByTagName('dl:documentId')
                    source = str(source[0].firstChild.data)
                    reactionSmarts = reaction.getElementsByTagName('dl:reactionSmiles')
                    for reactionSmart in reactionSmarts:
                        reactionSmart = str(reactionSmart.firstChild.data.split(' ')[0])
                        r,s,p = reactionSmart.split('>')
                        r = r.split('.')
                        p = p.split('.')
                        for i in r:
                            if 'C:' in i or 'CH:' in i or 'CH2:' in i or 'CH3:' in i:
                                rlist.append(i)
                            else:
                                if ":" in i:
                                    try: 
                                        i = remove_mapping(i)
                                    except Exception as e:
                                        print(e)
                                        i = i
                                else:
                                    try: 
                                        i = Chem.CanonSmiles(i)
                                    except Exception as e:
                                        print(e)
                                        i = i
                                rglist.append(i)
                    rlist,rglist = check_charge(rlist,rglist)
                    products = reaction.getElementsByTagName("product")
                    plist = []
                    for product in products:
                        identifiers = product.getElementsByTagName("identifier")
                        for identifier in identifiers:
                            if identifier.getAttribute('dictRef') == "cml:smiles":
                                #print(identifier.getAttribute('value'))
                                plist.append(identifier.getAttribute('value'))
                     

                    reactants = reaction.getElementsByTagName("reactant")
                    all_rlist = [] # USPTO considers reactants and reagents as reactants,so we need to split them.
                    all_rglist = []
                    for reactant in reactants:
                        identifiers = reactant.getElementsByTagName("identifier")
                        for identifier in identifiers:
                            if identifier.getAttribute('dictRef') == "cml:smiles":
                                #print(identifier.getAttribute('value'))
                                mol = identifier.getAttribute('value')
                                mols = mol.split('.')
                                isrg = True
                                for i in mols:
                                    try: 
                                        i = Chem.CanonSmiles(i)
                                    except Exception as e:
                                        print(e)
                                        i = i
                                    if i not in rglist:
                                        isrg = False
                                        break
                                if isrg:
                                    all_rglist.append(mol)
                                else:
                                    all_rlist.append(mol)
                    if len(all_rglist) == 0:
                        all_rglist = 'None'
                    


                    spectators = reaction.getElementsByTagName("spectator")
                    slist = [] # solvent list
                    clist = [] # catalyst list
                    for spectator in spectators:
                        identifiers = spectator.getElementsByTagName("identifier")
                        for identifier in identifiers:
                            if identifier.getAttribute('dictRef') == "cml:smiles":
                                #print(identifier.getAttribute('value'))
                                if spectator.getAttribute('role') == "solvent":
                                    slist.append(identifier.getAttribute('value'))
                                else:
                                    clist.append(identifier.getAttribute('value'))

                    if len(plist) == 0 or len(all_rlist) == 0:
                        continue
                    alist.append(n)
                    alist.append(source)
                    alist.append(reactionSmart)
                    alist.append(plist)
                    alist.append(all_rlist)
                    alist.append(all_rglist)
                    if len(clist) == 0:
                        alist.append('None')
                    else:
                        alist.append(clist)
                    if len(slist) == 0:
                        alist.append('None')
                    else:
                        alist.append(slist)
                    
                    
                    alllist.append(alist)
                    n +=1

            writer.writerows(alllist)
            print('-------------------------------------------')
            print('%s done'%filename)
            print('-------------------------------------------')
    all_data = pd.read_csv("%s/%s_all.csv"%(args.save_path,args.data_name))
    all_data.drop_duplicates(inplace=True,subset =['reaction','products','reactants','reagent','catalyst','solvent'])
    all_data.to_csv(out_path)
    print('alldone')


def get_m_data(args):
    global data
    '''
    This function is used to get the data that the template has more than m reactions.
    data_name: the name of the dataset
    save_path: the path to save the results
    m: the minimum number of occurrences for a reaction type.
    '''
    data_path = "%s/%s_all.csv"%(args.save_path,args.data_name)
    temp_path = "%s/classif_by_temp.csv"%args.save_path
    out_path = "%s/%s_%s+.csv"%(args.save_path,args.data_name,args.min_num_covered_rxns_by_rxn_centralized_template)
    n = 0
    all_data = []
    data = pd.read_csv(data_path)
    with open(temp_path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            lens = len(classes)  
            print("num of temp:",lens)
            for tem in classes:
                if len(eval(classes[tem])) < args.min_num_covered_rxns_by_rxn_centralized_template:
                    continue
                if tem == 'else':
                    continue
                if n%(lens//1000) == 0:
                    print("\r", end="")
                    print("Progress: {}%: ".format(n/((lens//1000)*10)), "▋" * int(n/((lens//1000)*20)), end="")
                    sys.stdout.flush()
                go = False
                for i in eval(classes[tem]):
                    adic = {}
                    adic['_id'] = data['_id'][int(i)-1]
                    adic['source'] = data['source'][int(i)-1]
                    adic['reaction'] = data['reaction'][int(i)-1]
                    adic['products'] = '.'.join(i for i in eval(data['products'][int(i)-1]))
                    adic['reactants'] = '.'.join(i for i in eval(data['reactants'][int(i)-1]))
                    adic['template'] = n
                    adic['tem_smart'] = tem
                    cats = data['catalyst'][int(i)-1]
                    solvs = data['solvent'][int(i)-1]
                    reag = data['reagent'][int(i)-1]
                    if cats == 'None':
                        cat = cats
                    else:
                        cats = Cleam_merge(cats)
                        cat = '.'.join(i for i in cats)
                    adic['cat'] = cat
                    if solvs == 'None':
                        solvs = solvs
                    else:
                        solvs = Cleam_merge(solvs)
                    if solvs == 'None':
                        for i in range(2):
                            adic['solv%s'%i] = 'None'
                    else:
                        for i in range(2):
                            if len(solvs) > i:
                                adic['solv%s'%i] = solvs[i]
                            else:
                                adic['solv%s'%i] = 'None'
                    if reag == 'None':
                        for i in range(4):
                            adic['reag%s'%i] = 'None'
                    else:
                        reag = Cleam_merge(reag)
                        for i in range(4):
                            if len(reag) > i:
                                adic['reag%s'%i] = reag[i]
                            else:
                                adic['reag%s'%i] = 'None'
                    go = True
                    all_data .append(adic)
                if go:
                    n +=1
                    
                    
    all_data = pd.DataFrame(all_data)
    all_data.drop_duplicates(inplace=True)
    all_data.to_csv(out_path)
    print('done')

def split_data(args):
    cat_list_N,solv_list_N,reag_list_N = open_csv(args.save_path)
    data_path = "%s/%s_%s+.csv"%(args.save_path,args.data_name,args.min_num_covered_rxns_by_rxn_centralized_template)
    condition_list_path = "%s/condition_list.json.gz"%args.save_path
    with gzip.open(condition_list_path) as f:
        condition_list = json.load(f)
    data = pd.read_csv(data_path)
    drop_list = []
    for i in range(data.shape[0]):
        tem = data['template'][i]
        conditions = [data['cat'][i],data['solv0'][i],data['solv1'][i],data['reag0'][i],data['reag1'][i],data['reag2'][i]]
        for j in range(len(conditions)):
            if str(conditions[j]) == 'nan':
                conditions[j] = "None"
        if conditions[0] not in cat_list_N:
            drop_list.append(False)
        elif conditions[1] not in solv_list_N or conditions[2] not in solv_list_N:
            drop_list.append(False)
        elif conditions[3] not in reag_list_N or conditions[4] not in reag_list_N or conditions[5] not in reag_list_N:
            drop_list.append(False)
        elif condition_list[tem]['conditions'][str(conditions)] < 2:
            drop_list.append(False)
        else:
            drop_list.append(True)
    data = data[drop_list]
    data_train = data.sample(frac=0.90,random_state=1)
    data_test_val = data.drop(data_train.index)
    data_test = data_test_val.sample(frac=0.5,random_state=1)
    data_val = data_test_val.drop(data_test.index)
    data_train.to_csv("%s/data_train.csv"%(args.save_path),index=False)
    data_test.to_csv("%s/data_test.csv"%(args.save_path),index=False)
    data_val.to_csv("%s/data_val.csv"%(args.save_path),index=False)
    print('done')



if __name__ == '__main__':
    pass