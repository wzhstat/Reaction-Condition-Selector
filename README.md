# Reaction-Condition-Selector
This repository contains a reaction condition selector.

# Python Version
python 3.10

# Install Requirements
To run the reaction condition selector you need:
* RDkit
* RDchiral <br>

You can go to https://github.com/connorcoley/rdchiral for more informaton.

# Get Dataset
Our data is derived from a dataset extracted from patents provided by the USPTO, which can be downloaded from the following address: 
https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873 <br>
We need two of these sections, ```1976_Sep2016_USPTOgrants_cml.7z``` and ```2001_Sep2016_USPTOapplications_cml.7z```. After decompression, plaese save the file under applications directly to the grants file.<br>

## Config Settings
To get the data for training, you can run ```get_data.py```, this  uses only the CPU and takes about 6 hours to run on a laptop. 
```
python get_data.py --data_name USPTO --data_path ./data/grants --save_path ./data
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
```USPTO_all.csv```: It contains information extracted directly from the USPTO data set, includes *id*,*source*,*reaction*, *products*, *reactants*, *reagent*, *catalyst*, *solvent*. Reaction is stored in Smart format, products, reactants, reagent, catalyst and solvent are stored in list form.<br>

```
, _id, source, reaction, products, reactants, reagent, catalyst, solvent
0, 1, US03930836, [Br:1][CH2:2][CH2:3][OH:4].[CH2:5]([S:7](Cl)(=[O:9])=[O:8])[CH3:6].CCOCC>C(N(CC)CC)C>[CH2:5]([S:7]([O:4][CH2:3][CH2:2][Br:1])(=[O:9])=[O:8])[CH3:6], ['C(C)S(=O)(=O)OCCBr'], ['BrCCO', 'C(C)S(=O)(=O)Cl']	['CCOCC'], None, ['C(C)N(CC)CC']
1, 2, US03930836, [Br:1][CH2:2][CH2:3][CH2:4][OH:5].[CH3:6][S:7](Cl)(=[O:9])=[O:8].CCOCC>C(N(CC)CC)C>[CH3:6][S:7]([O:5][CH2:4][CH2:3][CH2:2][Br:1])(=[O:9])=[O:8], ['CS(=O)(=O)OCCCBr'], ['BrCCCO', 'CS(=O)(=O)Cl'], ['CCOCC'], None, ['C(C)N(CC)CC']
2, 3, US03930836, [CH2:1]([Cl:4])[CH2:2][OH:3].CCOCC.[CH2:10]([S:14](Cl)(=[O:16])=[O:15])[CH:11]([CH3:13])[CH3:12]>C(N(CC)CC)C>[CH2:10]([S:14]([O:3][CH2:2][CH2:1][Cl:4])(=[O:16])=[O:15])[CH:11]([CH3:13])[CH3:12], ['C(C(C)C)S(=O)(=O)OCCCl'], ['C(CO)Cl', 'C(C(C)C)S(=O)(=O)Cl'], ['CCOCC'], None, ['C(C)N(CC)CC']
6, 7, US03930837, [Cl:1][C:2]1[N:3]=[CH:4][C:5]2[C:10]([CH:11]=1)=[C:9]([N+:12]([O-])=O)[CH:8]=[CH:7][CH:6]=2.O.[OH-].[Na+]>C(O)(=O)C.[Fe]>[Cl:1][C:2]1[N:3]=[CH:4][C:5]2[C:10]([CH:11]=1)=[C:9]([NH2:12])[CH:8]=[CH:7][CH:6]=2, ['ClC=1N=CC2=CC=CC(=C2C1)N'], ['ClC=1N=CC2=CC=CC(=C2C1)[N+](=O)[O-]'], ['O', '[OH-].[Na+]'], ['[Fe]'], ['C(C)(=O)O']
```

```USPTO.csv```: It is the data obtained after weight removal.<br>
```USPTO_m+.csv```: It contains the data that the template has more than m reactions.<br>
```classif_by_temp```: The file shows the corresponding reaction for each template.<br>
```all_cat_withoutN.csv```, ```all_solv_withN.csv```, ```all_solv_withoutN.csv```, ```all_solv_withN.csv```, ```all_reag_withoutN.csv```, ```all_reag_withN.csv```: These are statistics on reaction conditions.<br>
```data_train.csv```,```data_test.csv```,```data_val.csv```:These are the data sets that are randomly divided.<br>
The processed data can also be obtained from https://www.dropbox.com/scl/fo/v1rhyes2wvead9dz3x4fb/h?rlkey=nqtst7azldcry3ixnoigmcv3v&dl=0 <br>

# Training D-MPNN Model

## Step 1 Data Pre-processing
Before training, make sure you have the data file, which should contain two parts: the csv files```data_train.csv```, ```data_test.csv``` and ``` data_val.csv``` and the ```keys``` folder.<br>
You can run the ```preprocessing.py``` file to get the preprocessed data.<br>
```
python preprocessing.py
```


## Step 2 Training
You can run ```train.sh``` files directly to get models of catalyst, solvents, and reagents, we recommend using GPUs for faster training. Corresponding models are also given in Models. <br>
To train model for a particular condition, take solv0 for example, you can run:<br>
```
chemprop_train --target_columns solv0 --data_path ./data/MPNN_data/GCN_data/GCN_data_train.csv  --separate_val_path ./data/MPNN_data/GCN_data/GCN_data_val.csv --separate_test_path ./data/MPNN_data/GCN_data/GCN_data_test.csv --dataset_type multiclass --multiclass_num_classes 538 --save_dir ./data/models/GCN_solv0  --reaction --extra_metrics accuracy top3 --epochs 35
```
<br>
You end up with a models folder that contains the trained D-MPNN model. You can also download trained models directly from models/models. <br>


# Using Trained Model to Make Predictions

## Step 1 Construction of Reaction Condition Library
Before making a prediction, run ```BuildConditionLibrary.sh``` to build the three types of reaction condition libraries needed for prediction, namely the r1 library, the r0 library, and the r0* library. You can also run get_condition_library.py directly to get a specific library of reaction conditions.
```
python class_conditions.py --Inclusion 0.8 --data_set train --tpl_radius 0
```
```Inclusion```is the tolerance for labels of each category.   It indicates that a label can be selected as a category label only when the number of times it appears is greater than the total number of conditions times Inclusion.<br>
```data_set``` is the data set to be collected. <br>
```--tpl_radius``` is the radius of templates used for categorization.<br>

## Step 2 Prediction
You can run ```make_predictions.sh``` to complete the model predictions, which generates a prediction folder containing predictions with and without clustering.<br>
An examples of prediction results that include clustering are as follows:<br>
```
class id: 50_4
class label [[], ['phosphine', 'halide']]
best condition: ['None', 'None', 'None', 'BrP(Br)Br', 'None', 'None']
score: 0.43502089128291327
condition score: [["['None', 'None', 'None', 'BrP(Br)Br', 'None', 'None']", 0.43502089128291327], ["['None', 'None', 'None', 'O=P(Br)(Br)Br', 'None', 'None']", 0.4005131881456007], ["['None', 'CN(C)C=O', 'None', 'BrP(Br)Br', 'None', 'None']", 0.02543110079025514], ["['None', 'CN(C)C=O', 'None', 'O=P(Br)(Br)Br', 'None', 'None']", 0.023413797956965517]... ["['None', 'None', 'None', 'CCOC(C)=O.O=C([O-])O.[Na+]', 'O=P(Cl)(Cl)Cl', 'None']", 6.239001680815093e-37]]
--------------------------------------------
class id: 50_5
class label [[], ['halide', 'aldehyde', 'amine']]
best condition: ['None', 'CN(C)C=O', 'None', 'O=S(Cl)Cl', 'None', 'None']
score: 6.458456286464104e-14
condition score: [["['None', 'CN(C)C=O', 'None', 'O=S(Cl)Cl', 'None', 'None']", 6.458456286464104e-14], ["['None', 'None', 'None', 'CN(C)C=O', 'O=S(Cl)Cl', 'None']", 1.2491413029613758e-18], ["['None', 'ClCCl', 'None', 'CN(C)C=O', 'O=S(Cl)Cl', 'None']", 7.428030276709379e-22], ["['None', 'Cc1ccccc1', 'None', 'CN(C)C=O', 'O=S(Cl)Cl', 'None']", 3.421287350994334e-22], ["['None', 'ClCCl', 'None', 'CN(C)C=O', 'O=C(Cl)C(=O)Cl', 'None']", 1.1259136045695308e-23],...["['None', 'CN(C)C=O', 'None', 'O=C(Cl)C(=O)Cl', 'None', 'None']", 7.935796084468542e-26], ["['None', 'CN(C)C=O', 'None', 'CCN(C=O)CC', 'O=S(Cl)Cl', 'None']", 4.955859832749846e-32]]

...

--------------------------------------------
class id: 50_0
class label [['ionic'], ['phosphine', 'halide']]
best condition: ['C[N+](C)(C)C.[Cl-]', 'None', 'None', 'O=P(Cl)(Cl)Cl', 'None', 'None']
score: 9.45609355368023e-28
condition score: [["['C[N+](C)(C)C.[Cl-]', 'None', 'None', 'O=P(Cl)(Cl)Cl', 'None', 'None']", 9.45609355368023e-28], ["['CC[N+](CC)(CC)CC.[Cl-]', 'CC#N', 'None', 'CN(C)c1ccccc1', 'O=P(Cl)(Cl)Cl', 'None']", 4.047334473315928e-40]]
```
## Setp 3 Calculate Accuracy
Please run Score.py to calculate the accuracy of the model predictions.
```
python Score.py --data_path ./data/data_test.csv --pred_path ./data/prediction
```

# Predictions of conditions in actual drug synthesis routes
## Step 1 Obtaining route data
Our actual route data was artificially extracted from the article, and the processed data file can be downloaded from https://www.dropbox.com/scl/fo/v1rhyes2wvead9dz3x4fb/h?rlkey=nqtst7azldcry3ixnoigmcv3v&dl=0.<Br>


