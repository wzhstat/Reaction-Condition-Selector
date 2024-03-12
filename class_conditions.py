from ConditionClassifier import get_temp_condition, Classify_reaction_conditions
import argparse
from joblib import Parallel, delayed
import json

def main():
    parser = argparse.ArgumentParser(description='Classify reaction conditions')
    parser.add_argument('--data_path', type=str, default='data', help='path to data')
    parser.add_argument('--data_name', type=str, default='USPTO', help='name of data')
    parser.add_argument('--save_path', type=str, default='data', help='path to save')
    parser.add_argument('--min_num_covered_rxns_by_rxn_centralized_template', type=int, default=5)
    parser.add_argument('--Inclusion', type=int, default=0.8, help='Inclusion of reaction conditions label')
    parser.add_argument('--data_set', type=str, default='train', help='data set to use,can be all, train, test or val.')
    parser.add_argument('--n_jobs', type=int, default=-1, help='number of jobs to run in parallel')
    parser.add_argument('--tpl_radius', type=int, default=0, help='radius of template')
    args = parser.parse_args()
    condition_list = get_temp_condition(args)
    print('start to classify conditions')
    classed_condition_list = Parallel(n_jobs=-1, verbose=4)(delayed(Classify_reaction_conditions)(sorted(condition_list[i]['conditions']),condition_list[i]['tpl'],condition_list[i]['tpl_smarts'],args) for i in range(len(condition_list)))
    classed_condition_dic = {}
    for i in range(len(classed_condition_list)):
        classed_condition_dic[condition_list[i]['tpl_smarts']] = classed_condition_list[i]
    with open('%s/classed_conditions_library.json'%(args.save_path),'w') as f:
        json.dump(classed_condition_dic,f)
    print('classify conditions done')

if __name__ == '__main__':
    main()
