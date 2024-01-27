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
To get the data for training, you can run ```get_data.py```, this  uses only the CPU and takes about 6 hours to run on a laptop. 
```
python get_data.py --data_name USPTO --data_path datapath --save_path ./data/
```
Before runing the program, here are some configs that you can change:
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

```
	_id	source	reaction	products	reactants	reagent	catalyst	solvent
0	1	US03930836	[Br:1][CH2:2][CH2:3][OH:4].[CH2:5]([S:7](Cl)(=[O:9])=[O:8])[CH3:6].CCOCC>C(N(CC)CC)C>[CH2:5]([S:7]([O:4][CH2:3][CH2:2][Br:1])(=[O:9])=[O:8])[CH3:6]	['C(C)S(=O)(=O)OCCBr']	['BrCCO', 'C(C)S(=O)(=O)Cl']	['CCOCC']	None	['C(C)N(CC)CC']
1	2	US03930836	[Br:1][CH2:2][CH2:3][CH2:4][OH:5].[CH3:6][S:7](Cl)(=[O:9])=[O:8].CCOCC>C(N(CC)CC)C>[CH3:6][S:7]([O:5][CH2:4][CH2:3][CH2:2][Br:1])(=[O:9])=[O:8]	['CS(=O)(=O)OCCCBr']	['BrCCCO', 'CS(=O)(=O)Cl']	['CCOCC']	None	['C(C)N(CC)CC']
2	3	US03930836	[CH2:1]([Cl:4])[CH2:2][OH:3].CCOCC.[CH2:10]([S:14](Cl)(=[O:16])=[O:15])[CH:11]([CH3:13])[CH3:12]>C(N(CC)CC)C>[CH2:10]([S:14]([O:3][CH2:2][CH2:1][Cl:4])(=[O:16])=[O:15])[CH:11]([CH3:13])[CH3:12]	['C(C(C)C)S(=O)(=O)OCCCl']	['C(CO)Cl', 'C(C(C)C)S(=O)(=O)Cl']	['CCOCC']	None	['C(C)N(CC)CC']
6	7	US03930837	[Cl:1][C:2]1[N:3]=[CH:4][C:5]2[C:10]([CH:11]=1)=[C:9]([N+:12]([O-])=O)[CH:8]=[CH:7][CH:6]=2.O.[OH-].[Na+]>C(O)(=O)C.[Fe]>[Cl:1][C:2]1[N:3]=[CH:4][C:5]2[C:10]([CH:11]=1)=[C:9]([NH2:12])[CH:8]=[CH:7][CH:6]=2	['ClC=1N=CC2=CC=CC(=C2C1)N']	['ClC=1N=CC2=CC=CC(=C2C1)[N+](=O)[O-]']	['O', '[OH-].[Na+]']	['[Fe]', '[Fe]']	['C(C)(=O)O']

```





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

