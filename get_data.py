from DataProcessor import extract_data, count_common_data, get_tem
from condition_classifier import get_temp_condition
import time
import argparse
import pandas as pd


if __name__ == '__main__':
    import site
    print(site.getusersitepackages())
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_name', type=str, default='USPTO')
    parser.add_argument('--data_path', type=str, default='./data/grants')
    parser.add_argument('--save_path', type=str, default='./data')
    parser.add_argument('--min_num_covered_rxns_by_rxn_centralized_template', type=int, default=5)
    parser.add_argument('--min_num_covered_rxns_by_solvent', type=int, default=5)
    parser.add_argument('--min_num_covered_rxns_by_reagent', type=int, default=5)
    parser.add_argument('--min_num_covered_rxns_by_catalyst', type=int, default=5)
    args = parser.parse_args()


    t0 = time.time()
    print("------------------")
    print("start to get data")
    extract_data.get_datas(args)
    t1 = time.time()
    print("extract data done")
    print("time used:", t1-t0)
    print("------------------")
    print("start to get template")
    get_tem.classif_by_temp(args)
    t2 = time.time()
    print("get template done")
    print("time used:", t2-t1)
    print("------------------")
    print("start to statistic data")
    count_common_data.count_data(args)
    t3 = time.time()
    print("statistic data done")
    print("time used:", t3-t2)
    print("------------------")
    print("start to get m data")
    extract_data.get_m_data(args)
    t4 = time.time()
    print("get m data done")
    print("time used:", t4-t3)
    print("------------------")
    print("start to get condition")
    get_temp_condition(args)
    t5 = time.time()
    print("get condition done")
    print("time used:", t5-t4)
    print("------------------")
    print("start to split data")
    extract_data.split_data(args)
    t6 = time.time()
    print("split data done")
    print("time used:", t6-t5)