from data_processor import ExtractData, CountCommonData, GetTem

config = {
    'data_name':'1976-2016',
    'data_path':'./data/grants',
    'save_path':'./data',
    "min_num_covered_rxns_by_rxn_centralized_template":5,
    "min_num_covered_rxns_by_solvent":5,
    "min_num_covered_rxns_by_reagent":5,
    "min_num_covered_rxns_by_catalyst":5,

}



if __name__ == '__main__':
    print("------------------")
    print("start to get data")
    extract_data.get_datas(config['data_name'], config['data_path'], config['save_path'])
    print("extract data done")
    print("------------------")
    print("start to get template")
    get_tem.classif_by_temp(config['data_name'], config['save_path'])
    print("get template done")
    print("------------------")
    print("start to get m data")
    extract_data.get_m_data(config['data_name'], config['save_path'], config['min_num_covered_rxns_by_rxn_centralized_template'])
    print("get m data done")
    print("------------------")
    print("start to statistic data")
    count_common_data.count_data(config['data_name'], config['save_path'], config['min_num_covered_rxns_by_rxn_centralized_template'], config['min_num_covered_rxns_by_catalyst'], config['min_num_covered_rxns_by_solvent'], config['min_num_covered_rxns_by_reagent'])
    print("statistic data done")
    print("------------------")



