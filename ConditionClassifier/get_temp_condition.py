import pandas as pd
import numpy as np
import sys


def get_temp_condition(args):
    '''
    This function is used to count the reaction conditions contained under each template.
    '''
    print('start to get condition list')
    if args.data_set == 'all':
        data =pd.read_csv('%s/%s_%s+.csv'%(args.data_path,args.data_name,args.min_num_covered_rxns_by_rxn_centralized_template))
    elif args.data_set == 'train':
        data =pd.read_csv('%s/data_train.csv'%(args.data_path))
    elif args.data_set == 'test':
        data =pd.read_csv('%s/data_test.csv'%(args.data_path))
    elif args.data_set == 'val':
        data =pd.read_csv('%s/data_val.csv'%(args.data_path))
    data.sort_values(by=['template_r%s'%args.tpl_radius],inplace=True)
    l = data.shape[0]
    condition_list = {}
    for i in range(l):
        if i%(l//1000) == 0:
            print("\r", end="")
            print("Progress: {}%: ".format(i/((l//1000)*10)), "▋" * int(i/((l//1000)*20)), end="")
            sys.stdout.flush()
        tem = str(data['template_r%s'%args.tpl_radius][i])
        tem_smarts = data['tpl_SMARTS_r%s'%args.tpl_radius][i]
        cat = data['cat'][i] if str(data['cat'][i]) != 'nan' else 'None'
        solv0 = data['solv0'][i] if str(data['solv0'][i]) != 'nan' else 'None'
        solv1 = data['solv1'][i] if str(data['solv1'][i]) != 'nan' else 'None'
        reag0 = data['reag0'][i] if str(data['reag0'][i]) != 'nan' else 'None'
        reag1 = data['reag1'][i] if str(data['reag1'][i]) != 'nan' else 'None'
        reag2 = data['reag2'][i] if str(data['reag2'][i]) != 'nan' else 'None'
        reactants = data['reactants'][i]
        products = data['products'][i]
        if tem_smarts in condition_list:
            condition = [cat,solv0,solv1,reag0,reag1,reag2]
            if str(condition) in condition_list[tem_smarts]['conditions']:
                condition_list[tem_smarts]['conditions'][str(condition)].append((reactants,products))
            else:
                condition_list[tem_smarts]['conditions'][str(condition)] = [(reactants,products)]
        else:
            condition_list[tem_smarts] = {'tpl':str(tem),'tpl_smarts':data['tpl_SMARTS_r%s'%args.tpl_radius][i],'conditions':{}}
            condition = [cat,solv0,solv1,reag0,reag1,reag2]
            if str(condition) in condition_list[tem_smarts]['conditions']:
                condition_list[tem_smarts]['conditions'][str(condition)].append((reactants,products))
            else:
                condition_list[tem_smarts]['conditions'][str(condition)] = [(reactants,products)]
    print('get condition list done')
    return condition_list

def default_dump(obj):
    """Convert numpy classes to JSON serializable objects."""
    if isinstance(obj, (np.integer, np.floating, np.bool_)):
        return obj.item()
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj

if __name__ == '__main__':
    #get_temp_condition()
    pass

