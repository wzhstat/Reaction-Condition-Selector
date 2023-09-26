import os
from xml.dom import minidom
from xml.dom.minidom import parse
import csv
import csv
import pandas as pd
from rdkit import Chem
from rdchiral import template_extractor

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

def Cleam_merge(alist):
    '''
    This function is used to clean the smiles and merge the same smiles.
    '''
    newlist = eval(str(alist))
    outlist = []
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
            outlist.append(s)
        outlist.sort()
    return outlist

def check_charge(rlist,rglist):
    '''
    This function is used to check the charge of the reactants and reagents.
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



def get_datas(load_path="./data/grants",out_path="./data/1976-20161.csv"):
    '''
    This function is used to get the data from the xml files and save them into a csv file.
    load_path: the path of the data
    save_path: the path of the csv file
    '''
    header = ['_id','reaction','products','reactants','reagent','catalyst','solvent']
    n = 1
    with open(out_path, 'w', newline='') as file:
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
                    
                    reactionSmarts = reaction.getElementsByTagName('dl:reactionSmiles')
                    for reactionSmart in reactionSmarts:
                        reactionSmart = str(reactionSmart.firstChild.data.split(' ')[0])
                        r,s,p = reactionSmart.split('>')
                        r = r.split('.')
                        p = p.split('.')
                        for i in r:
                            if ':' in i:
                                rlist.append(i)
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

                    
                    alist.append(n)
                    alist.append(reactionSmarts)
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
    print('alldone')


def get_5_data(data_path = "data/1976-2016.csv", temp_path = "data/classif_by_temp.csv", out_path = "data/data_5+.csv"):
    '''
    This function is used to get the data that the template has more than 5 reactions.
    '''
    n = 0
    all_data = []
    data = pd.read_csv(data_path)
    with open(temp_path,'r') as f:
        reader = csv.DictReader(f)
        for classes in reader:
            lens = len(eval(classes))
            for tem in classes:
                if n%(lens//100) == 0:
                    print(n/(lens//100),'%')
                go = False
                if len(eval(classes[tem])) < 5:
                    continue
                for i in eval(classes[tem]):
                    adic = {}
                    cats = data['catalyst'][int(i)-1]
                    solvs = data['solvent'][int(i)-1]
                    reag = data['reagent'][int(i)-1]
                    if cats == 'None':
                        cat = cats
                    else:
                        cats = Cleam_merge(cats)
                        cat = '.'.join(i for i in cats)
                    if solvs == 'None':
                        solv = solvs
                    else:
                        solvs = Cleam_merge(solvs)
                        solv = '.'.join(i for i in solvs)
                    adic['reaction'] = data['reaction'][int(i)-1]
                    adic['products'] = '.'.join(i for i in eval(data['products'][int(i)-1]))
                    adic['reactants'] = '.'.join(i for i in eval(data['reactants'][int(i)-1]))
                    adic['template'] = n
                    adic['tem_smart'] = tem
                    adic['cat'] = cat
                    adic['solv'] = solv
                    if reag == 'None':
                        pass
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
    print(all_data)
    print(len(all_data))
    all_data.to_csv(out_path)
    print('done')
