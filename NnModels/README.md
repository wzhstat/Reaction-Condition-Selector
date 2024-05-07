# nn Models
This section contains some simple feedforward neural network models with inputs of reactants (512), product fingerprints (512) and other fingerprints(512*n) and outputs of vectors directed to the target.

## Models
```nnModel0```: An MLP model with two hidden layers, by default, n1 is 128 and n2 is 32. <br>
```nnModel1```: An MLP model with one full connection layer and two highway layers, n1 is 128 by default.<br>
```nnModel2```: An MLP model with one full connection layer and two highway layers, n1 is 128 by default, This model differs from Model 1 in that its input contains a reaction template in the form of one-hot in addition to the molecular fingerprints of the reactants and products<br>

## Training 
To train a model, you can run ```train_nnmodel.py```. Depending on the number of epochs selected, the training time is about 1-20h, but plaase note that using nnModel2 or ```condition filter``` will greatly increase running time. Here are some configs that you can change: <pr>
```
config = {
    'data_name':'1976-2016_5+',
    'save_path':'./data',
    'model_path':'./models',
    'model': MLPModel.nnModel1,
    'input': 'rfpgen+pfpgen+rxnfp', #['rfpgen+pfpgen','rfpgen+pfpgen+rxnfp','rfpgen+pfpgen+tem'] 
    'target':['cat','solv','reag0','reag1'],
    'withN': False,
    'epochs': { 'cat': 40, 'solv':10 , 'reag0': 30, 'reag1': 30},
    'n1': 128,
    'n2': 64,
    'Ir': 0.0001,
    'batch_size': 128,
    'condition filter': True,
    'Hierarchical prediction':True
}
```
* ```withN```: Indicates whether to filter None data from data
* ```input```: This variable is used to determine the input of the reaction, the conventional inputs are 'rfpgen+pfpgen','rfpgen+pfpgen+rxnfp','rfpgen+pfpgen+rxnfp','rfpgen+pfpgen+ pfpgen+tem','rfpgen+pfpgen+tem','rfpgen+pfpgen+ tem ', You can also control whether conditional information is added by adding '+cat', '+ solv ', '+reag0', etc. to the tail.
* ```target```: A list of the models you want to train.
* ```n1```: The size of the first hidden layer.
* ```n2```: The size of the second hidden layer.
* ```condition filter```: This parameter controls whether to add a result filter, which first counts the condition used by all templates, outputs a template condition library, and filters out the output that is not included in the library.
* ```Hierarchical prediction```: This variable initiates a hierarchical model training, in which the first trained model variable is added to the next model input. In the example of config, the catalyst information is converted into a 512-dimensional fingerprint input when training the solv model and the catalyst information and solvent information input when training the reag0 model.


# Condition Classifier


