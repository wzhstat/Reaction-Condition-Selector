from DataProcessor import ExtractData, CountCommonData, GetTem
import time

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
    t0 = time.time()
    print("------------------")
    print("start to get data")
    ExtractData.get_datas(config['data_name'], config['data_path'], config['save_path'])
    t1 = time.time()
    print("extract data done")
    print("time used:", t1-t0)
    print("------------------")
    print("start to get template")
    GetTem.classif_by_temp(config['data_name'], config['save_path'])
    t2 = time.time()
    print("get template done")
    print("time used:", t2-t1)
    print("------------------")
    print("start to get m data")
    ExtractData.get_m_data(config['data_name'], config['save_path'], config['min_num_covered_rxns_by_rxn_centralized_template'])
    t3 = time.time()
    print("get m data done")
    print("time used:", t3-t2)
    print("------------------")
    print("start to statistic data")
    CountCommonData.count_data(config['data_name'], config['save_path'], config['min_num_covered_rxns_by_rxn_centralized_template'], config['min_num_covered_rxns_by_catalyst'], config['min_num_covered_rxns_by_solvent'], config['min_num_covered_rxns_by_reagent'])
    t4 = time.time()
    print("statistic data done")
    print("time used:", t4-t3)
    print("------------------")
    t5 = time.time()
    print("all done")
    print("time used:", t5-t0)



