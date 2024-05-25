from GetMpnnData import get_MPNN_data, save_csv
import argparse
import pandas as pd
import pandas as pd


def main():
    """
    This function is the main entry point for the script. It sets up command line arguments for data path, key path, save path, 
    a boolean flag 'N', data name, and target. It then iterates over a list of data types (validation, test, training) and a list 
    of targets (cat, solv0, solv1, reag0, reag1, reag2). For each combination of data type and target, it parses the command line 
    arguments and retrieves the corresponding data using the 'get_MPNN_data' function. If the target is 'cat', it also retrieves 
    the 'reaction' data. The retrieved data is then stored in a pandas DataFrame 'out'. Finally, it saves the DataFrame 'out' to a 
    CSV file using the 'save_csv' function.
    """
    parser = argparse.ArgumentParser(description='Get data for GCN')
    parser.add_argument('--data_path', type=str, default='data', help='path to data')
    parser.add_argument('--key_path', type=str, default='data/labels', help='path to data')
    parser.add_argument('--save_path', type=str, default='data/MPNN_data', help='path to save')
    parser.add_argument('--N', type=bool, default=True, help='whether to add None')
    parser.add_argument('--data_name', type=str)
    parser.add_argument('--target', type=str)
    data_list = ['data_val','data_test','data_train']
    target_list = ['cat','solv0','solv1','reag0','reag1','reag2']
    for data in data_list:
        out = pd.DataFrame()
        for target in target_list:
            args = parser.parse_args(['--data_name', data, '--target', target])
            if target == 'cat':
                out['reaction'] = get_MPNN_data(args)['reaction']
                out[target] = get_MPNN_data(args)['target']
            else:
                out[target] = get_MPNN_data(args)['target']
        save_csv(args,out)

if __name__ == '__main__':
    main()
