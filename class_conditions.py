from condition_classifier import get_temp_condition, Classify_reaction_conditions
import argparse
import json

def main():
    parser = argparse.ArgumentParser(description='Classify reaction conditions')
    parser.add_argument('--data_path', type=str, default='data', help='path to data')
    parser.add_argument('--data_name', type=str, default='USPTO', help='name of data')
    parser.add_argument('--save_path', type=str, default='data', help='path to save')
    parser.add_argument('--min_num_covered_rxns_by_rxn_centralized_template', type=int, default=5)
    parser.add_argument('--Inclusion', type=int, default=0.9, help='Inclusion of reaction conditions label')
    parser.add_argument('--data_set', type=str, default='all', help='data set to use,can be all or train.')
    args = parser.parse_args()
    condition_list = get_temp_condition(args)
    print('start to classify conditions')
    classed_condition_list = [Classify_reaction_conditions(sorted(condition_list['conditions'][i]),condition_list['tem'][i],condition_list['tem_smart'][i],args) for i in range(len(condition_list))]
    with open('%s/classed_condition_list.json'%args.save_path,'w') as f:
        json.dump(classed_condition_list,f)
    print('classify conditions done')
if __name__ == '__main__':
    main()