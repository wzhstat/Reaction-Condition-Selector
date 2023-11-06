# Reaction-Condition-Selector
This repository contains a reaction condition selector.

# Python Version
python 3.10

# Install Requirements
To run the reaction condition selector you need:
* RDkit
* RDchiral <br>

You can go to https://github.com/connorcoley/rdchiral for more informaton.

# Data Processor
This file is used to process the reaction condition data of USPTO1976-2016, and requires the 1976_Sep2016_USPTOgrants_cm file, which can be downloaded from the following address: https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873 <br>

## Config Settings
To get the data for training, you can run ```get_data.py```, this  uses only the CPU and takes about 3 hours to run on a laptop. Before runing the program, here are some configs that you can change:
```
config = {
    'data_name':'1976-2016',
    'data_path':'./data/grants',
    'save_path':'./data',
    "min_num_covered_rxns_by_rxn_centralized_template":5,
    "min_num_covered_rxns_by_solvent":5,
    "min_num_covered_rxns_by_reagent":5,
    "min_num_covered_rxns_by_catalyst":5,
}
```
* ```data_name```: Is the name of the data.
* ```data_path```: Is where you save the 1976_Sep2016_USPTOgrants_cm file.
* ```save_path```: Is where you want to save the files.
* ```min_num_covered_rxns_by_rxn_centralized_template```: Just like it's name it's the number of reactions that the template contains at least, templates containing fewer reactions than this value will be removed
* ```min_num_covered_rxns_by_solvent```: This is the minimum number of times a certain solvent appears
* ```min_num_covered_rxns_by_reagent```: This is the minimum number of times a certain reagent appears
* ```min_num_covered_rxns_by_catalyst```: This is the minimum number of times a certain catalyst appears<br>

## Output Files

The following files are generated when you run the program: <br>
```1976-2016.csv```: It contains information extracted directly from the USPTO data set, includes *id*, *reaction*, *products*, *reactants*, *reagent*, *catalyst*, *solvent*. Reaction is stored in Smart format, products, reactants, reagent, catalyst and solvent are stored in list form.<br>

    _id	reaction	products	reactants	reagent	catalyst	solvent
    1	[Br:1][CH2:2][CH2:3][OH:4].[CH2:5]([S:7](Cl)(=[O:9])=[O:8])[CH3:6].CCOCC>C(N(CC)CC)C>[CH2:5]([S:7]([O:4][CH2:3][CH2:2][Br:1])(=[O:9])=[O:8])[CH3:6]	['C(C)S(=O)(=O)OCCBr']	['BrCCO', 'C(C)S(=O)(=O)Cl']	['CCOCC']	None	['C(C)N(CC)CC']
    2	[N:1](=[C:3]([CH2:7][CH2:8][CH2:9][CH2:10][CH2:11][CH3:12])[C:4](=O)[CH3:5])O.[CH3:13][C:14](=O)[CH2:15][CH2:16][CH2:17][CH2:18][CH2:19][CH2:20][CH3:21].[N:23](Cl)=O.N=C=O>[Zn].C(O)(=O)C>[CH3:13][C:14]1[C:15]([CH2:16][CH2:17][CH2:18][CH2:19][CH2:20][CH3:21])=[N:23][C:4]([CH3:5])=[C:3]([CH2:7][CH2:8][CH2:9][CH2:10][CH2:11][CH3:12])[N:1]=1	['CC1=NC(=C(N=C1CCCCCC)C)CCCCCC']	['N(O)=C(C(C)=O)CCCCCC', 'CC(CCCCCCC)=O', 'N(=O)Cl']	['N=C=O']	['[Zn]']	['C(C)(=O)O']
    3	[CH3:1][O:2][C:3](=[O:17])[C:4]1[C:9]([N+:10]([O-])=O)=[CH:8][C:7]2[O:13][CH2:14][CH2:15][O:16][C:6]=2[CH:5]=1>[Pd].C(O)(=O)C>[CH3:1][O:2][C:3](=[O:17])[C:4]1[C:9]([NH2:10])=[CH:8][C:7]2[O:13][CH2:14][CH2:15][O:16][C:6]=2[CH:5]=1	['COC(C1=CC2=C(C=C1N)OCCO2)=O']	['COC(C1=CC2=C(C=C1[N+](=O)[O-])OCCO2)=O']	None	['[Pd]']	['C(C)(=O)O']




```1976-2016_m+.csv```: It contains the data that the template has more than m reactions.<br>
```classif_by_temp```: The file shows the corresponding reaction for each template.<br>
```all_cat_withoutN.csv```, ```all_solv_withN.csv```, ```all_solv_withoutN.csv```, ```all_solv_withN.csv```, ```all_reag_withoutN.csv```, ```all_reag_withN.csv```: These are statistics on reaction conditions.<br>

```classif_by_temp```and```all_cat_withoutN.csv```, ```all_solv_withN.csv```, ```all_solv_withoutN.csv```, ```all_solv_withN.csv```, ```all_reag_withoutN.csv```, ```all_reag_withN.csv``` can also be downloaded directly from data.zip.

# Popularity Baseline

This file contains two statistical models ```popularity_module0```, ```popularity_module1```,that output the most common conditions for all reactions and the most common conditions for a particular template.

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

# MPNN Models
This work uses chemprop's model, with only simple changes to the inputs and outputs in some parts. The changed chemprop can be download in ```https://github.com/wzhstat/Reaction-Condition-Selector/chemprop```.

## Data
The input data learned by the MPNN model has two columns, which are reaction Smarts and the grouping of the target conditions. You can get it by running ```get_mpnn_data.py``` ,Here are some configs that you can change: <br>
```
config = {
    "path": "./data",
    "file_name": '1976-2016_5+',
    "target": ['cat','solv','reag0'],
    "conditons": None,
    "condition_moduel": '',
    "N": True
}
```
Data samples: <br>
```
smarts	cat
[Cl:1][CH2:2][CH2:3][CH2:4][S:5](Cl)(=[O:7])=[O:6].[OH:16][CH2:17][C:18]([CH3:35])([CH3:34])[C@@H:19]([O:26][Si:27]([CH3:33])([CH3:32])[C:28]([CH3:31])([CH3:30])[CH3:29])/[CH:20]=[CH:21]/[C:22]([O:24][CH3:25])=[O:23]>>[Cl:1][CH2:2][CH2:3][CH2:4][S:5]([O:16][CH2:17][C:18]([CH3:35])([CH3:34])[C@@H:19]([O:26][Si:27]([CH3:33])([CH3:32])[C:28]([CH3:30])([CH3:29])[CH3:31])/[CH:20]=[CH:21]/[C:22]([O:24][CH3:25])=[O:23])(=[O:7])=[O:6]	2
Cl[S:2]([CH2:5][CH2:6][CH2:7][NH:8][C:9](=[O:11])[CH3:10])(=[O:4])=[O:3].[OH:19][CH2:20][C:21]([CH3:38])([CH3:37])[C@@H:22]([O:29][Si:30]([CH3:36])([CH3:35])[C:31]([CH3:34])([CH3:33])[CH3:32])/[CH:23]=[CH:24]/[C:25]([O:27][CH3:28])=[O:26]>>[C:9]([NH:8][CH2:7][CH2:6][CH2:5][S:2]([O:19][CH2:20][C:21]([CH3:38])([CH3:37])[C@@H:22]([O:29][Si:30]([CH3:36])([CH3:35])[C:31]([CH3:32])([CH3:34])[CH3:33])/[CH:23]=[CH:24]/[C:25]([O:27][CH3:28])=[O:26])(=[O:4])=[O:3])(=[O:11])[CH3:10]	2
[Cl:1][CH2:2][CH2:3][CH2:4][S:5](Cl)(=[O:7])=[O:6].[OH:16][CH2:17][C:18]([CH3:35])([CH3:34])[C@@H:19]([O:26][Si:27]([CH3:33])([CH3:32])[C:28]([CH3:31])([CH3:30])[CH3:29])/[CH:20]=[CH:21]/[C:22]([O:24][CH3:25])=[O:23]>>[Cl:1][CH2:2][CH2:3][CH2:4][S:5]([O:16][CH2:17][C:18]([CH3:35])([CH3:34])[C@@H:19]([O:26][Si:27]([CH3:33])([CH3:32])[C:28]([CH3:30])([CH3:29])[CH3:31])/[CH:20]=[CH:21]/[C:22]([O:24][CH3:25])=[O:23])(=[O:7])=[O:6]	2
Cl[S:2]([CH2:5][CH2:6][CH2:7][NH:8][C:9](=[O:11])[CH3:10])(=[O:4])=[O:3].[OH:19][CH2:20][C:21]([CH3:38])([CH3:37])[C@@H:22]([O:29][Si:30]([CH3:36])([CH3:35])[C:31]([CH3:34])([CH3:33])[CH3:32])/[CH:23]=[CH:24]/[C:25]([O:27][CH3:28])=[O:26]>>[C:9]([NH:8][CH2:7][CH2:6][CH2:5][S:2]([O:19][CH2:20][C:21]([CH3:38])([CH3:37])[C@@H:22]([O:29][Si:30]([CH3:36])([CH3:35])[C:31]([CH3:32])([CH3:34])[CH3:33])/[CH:23]=[CH:24]/[C:25]([O:27][CH3:28])=[O:26])(=[O:4])=[O:3])(=[O:11])[CH3:10]	2
```

## Training
You can run ```train.sh``` files directly to get models of cat, solv, and reag0, we recommend using GPUs for faster training. Corresponding models are also given in Models. <br>
To train model for a particular condition, take solv for example, you can run:<br>
```
chemprop_train --data_path ./data/GCN_solv_data_withN/GCN_data.csv --dataset_type multiclass --multiclass_num_classes 2973 --save_dir ./data/GCN_solv_withN_checkpoint --reaction --extra_metrics accuracy top3 --epochs 20
```

