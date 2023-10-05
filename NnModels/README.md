# NnModels
This section contains some simple feedforward neural network models with inputs of reactants (512) and product fingerprints (512) and outputs of vectors directed to the target.

## Training 
To train a model, you can run ```train_nnmodel.py```. Depending on the number of epochs selected, the training time is about 1h. Here are some config you can change: <pr>
```
config = {
    'data_name':'1976-2016_5+',
    'save_path':'./data',
    'model_path':'./models',
    'model': MLPModel.nnModel1,
    'target':['cat','solv','reag0','reag1'],
    'withN': False,
    'epochs': { 'cat': 30, 'solv':10 , 'reag0': 10, 'reag1': 30},
    'n1': 128,
    'n2': 32,
    'Ir': 0.0001,
    'batch_size': 128
}
```
