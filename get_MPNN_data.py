from GetMpnnData import get_MPNN_data
import argparse

def main():
    parser = argparse.ArgumentParser(description='Get data for GCN')
    parser.add_argument('--data_path', type=str, default='data', help='path to data')
    parser.add_argument('--save_path', type=str, default='data/MPNN_data', help='path to save')
    parser.add_argument('--data_name', type=str, default='data_train', help='name of data')
    parser.add_argument('--target', type=str, help='target name')
    parser.add_argument('--N', type=bool, default=True, help='whether to add None')
    
    data_list = ['data_train','data_test','data_val']
    target_list = ['cat','solv0','solv1','reag0','reag1','reag2']
    for data in data_list:
        for target in target_list:
            args = parser.parse_args(['--data_path', 'data', '--data_name', data, '--target', target])
            get_MPNN_data(args)

if __name__ == '__main__':
    main()