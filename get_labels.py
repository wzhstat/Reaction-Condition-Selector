import pandas as pd
import argparse
import os

def count_data(data):
    count = {}
    for mol in data:
        if mol in count:
            count[mol] += 1
        else:
            count[mol] = 1
    # Sort the dictionary by value
    count = dict(sorted(count.items(), key=lambda item: item[1], reverse=True))
    return count

def get_labels(args):
    data = pd.read_csv(args.data_path)
    catalysts = []
    for column in args.catalyst_column:
        catalysts += data[column].tolist()
    solvents = []
    for column in args.solvent_column:
        solvents += data[column].tolist()
    reagents = []
    for column in args.reagent_column:
        reagents += data[column].tolist()

    catalyst_count = count_data(catalysts)
    solvent_count = count_data(solvents)
    reagent_count = count_data(reagents)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(os.path.join(args.output_dir, 'cat_labels.csv'), 'w') as f:
        f.write('cat,count\n')
        for key, value in catalyst_count.items():
            f.write(f'{key},{value}\n')

    with open(os.path.join(args.output_dir, 'solv_labels.csv'), 'w') as f:
        f.write('sol,count\n')
        for key, value in solvent_count.items():
            f.write(f'{key},{value}\n')

    with open(os.path.join(args.output_dir, 'reag_labels.csv'), 'w') as f:
        f.write('reag,count\n')
        for key, value in reagent_count.items():
            f.write(f'{key},{value}\n')

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--data_path', type=str, default='test_template_extractor.csv')
    argparser.add_argument('--catalyst_column', type=list, default=['cat'])
    argparser.add_argument('--solvent_column', type=list, default=['solv0', 'solv1'])
    argparser.add_argument('--reagent_column', type=list, default=['reag0', 'reag1', 'reag2'])
    argparser.add_argument('--output_dir', type=str, default='labels')
    args = argparser.parse_args()  
    get_labels(args)





