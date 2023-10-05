from NnModels import MLPModel, TrainModel

config = {
    'data_name':'1976-2016_5+',
    'save_path':'./data',
    'model_path':'./models',
    'model': MLPModel.nnModel2,
    'target':['cat','solv','reag0','reag1'],
    'withN': False,
    'epochs': { 'cat': 1, 'solv':1 , 'reag0': 1, 'reag1': 1},
    'n1': 128,
    'n2': 32,
    'Ir': 0.0001,
    'batch_size': 128
}

if __name__ == '__main__':
    for target in config['target']:
        print('start to train %s model'%target)
        TrainModel.train(config['model'], 
                         config['save_path'], 
                         config['data_name'], 
                         config['withN'],
                         target, 
                         config['epochs'][target], 
                         config['n1'], 
                         config['n2'],
                         config['Ir'], 
                         config['batch_size'])  
        print('train %s model done'%target)
        print('------------------')
    print('all done')
