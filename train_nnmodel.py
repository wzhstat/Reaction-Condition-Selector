from NnModels import MLPModel, TrainModel

config = {
    'data_name':'1976-2016_5+',
    'save_path':'./data',
    'model_path':'./models',
    'model': MLPModel.nnModel1,
    'input': 'rfpgen+pfpgen+rxnfp', #['rfpgen+pfpgen','rfpgen+pfpgen+rxnfp','rfpgen+pfpgen+tem'] 
    'target':['cat','solv','reag0','reag1'],
    'withN': False,
    'epochs': { 'cat': 1, 'solv':2 , 'reag0': 2, 'reag1': 2},
    'n1': 128,
    'n2': 64,
    'Ir': 0.0001,
    'batch_size': 128,
    'condition filter': True,
    'Hierarchical prediction':True
}

if __name__ == '__main__':
    if config['condition filter']:
        print('start to sum the tem reaction conditions')
        TrainModel.sum_tem_condition(config['withN'],config['save_path'],config['data_name'])
        print('tem reaction conditions sum done')   
        pass
    for target in config['target']:
        print('start to train %s model'%target)
        TrainModel.train(config['input'],
                         config['model'], 
                         config['save_path'], 
                         config['data_name'], 
                         config['withN'],
                         target, 
                         config['epochs'][target], 
                         config['n1'], 
                         config['n2'],
                         config['Ir'], 
                         config['batch_size'],
                         config['Hierarchical prediction']
                         )
        print('train %s model done'%target)
        print('------------------')
        if config['Hierarchical prediction']:
            config['input'] += '+%s'%target
    print('all done')
