# Reaction-Condition-Selector
This repository contains a reaction condition selector.

# Python Version
python 3.10

# Install Requirements
To run the reaction condition selector you need:
* RDkit
* RDchiral <br>

You can go to https://github.com/connorcoley/rdchiral for more informaton.

# Data preprocessing
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

# Training scoring model
This work uses chemprop's model, with only simple changes to the inputs and outputs in some parts. The changed chemprop can be download in ```https://github.com/wzhstat/Reaction-Condition-Selector/chemprop```.

## Data
The input data learned by the MPNN model has two columns, which are reaction Smarts and the grouping of the target conditions. You can get it by running ```get_MPNN_data.py``` 
Data samples: <br>
```
reaction, cat
[Cl:1][CH2:2][CH2:3][CH2:4][S:5](Cl)(=[O:7])=[O:6].[OH:16][CH2:17][C:18]([CH3:35])([CH3:34])[C@@H:19]([O:26][Si:27]([CH3:33])([CH3:32])[C:28]([CH3:31])([CH3:30])[CH3:29])/[CH:20]=[CH:21]/[C:22]([O:24][CH3:25])=[O:23]>>[Cl:1][CH2:2][CH2:3][CH2:4][S:5]([O:16][CH2:17][C:18]([CH3:35])([CH3:34])[C@@H:19]([O:26][Si:27]([CH3:33])([CH3:32])[C:28]([CH3:30])([CH3:29])[CH3:31])/[CH:20]=[CH:21]/[C:22]([O:24][CH3:25])=[O:23])(=[O:7])=[O:6], 2
Cl[S:2]([CH2:5][CH2:6][CH2:7][NH:8][C:9](=[O:11])[CH3:10])(=[O:4])=[O:3].[OH:19][CH2:20][C:21]([CH3:38])([CH3:37])[C@@H:22]([O:29][Si:30]([CH3:36])([CH3:35])[C:31]([CH3:34])([CH3:33])[CH3:32])/[CH:23]=[CH:24]/[C:25]([O:27][CH3:28])=[O:26]>>[C:9]([NH:8][CH2:7][CH2:6][CH2:5][S:2]([O:19][CH2:20][C:21]([CH3:38])([CH3:37])[C@@H:22]([O:29][Si:30]([CH3:36])([CH3:35])[C:31]([CH3:32])([CH3:34])[CH3:33])/[CH:23]=[CH:24]/[C:25]([O:27][CH3:28])=[O:26])(=[O:4])=[O:3])(=[O:11])[CH3:10], 2
[Cl:1][CH2:2][CH2:3][CH2:4][S:5](Cl)(=[O:7])=[O:6].[OH:16][CH2:17][C:18]([CH3:35])([CH3:34])[C@@H:19]([O:26][Si:27]([CH3:33])([CH3:32])[C:28]([CH3:31])([CH3:30])[CH3:29])/[CH:20]=[CH:21]/[C:22]([O:24][CH3:25])=[O:23]>>[Cl:1][CH2:2][CH2:3][CH2:4][S:5]([O:16][CH2:17][C:18]([CH3:35])([CH3:34])[C@@H:19]([O:26][Si:27]([CH3:33])([CH3:32])[C:28]([CH3:30])([CH3:29])[CH3:31])/[CH:20]=[CH:21]/[C:22]([O:24][CH3:25])=[O:23])(=[O:7])=[O:6], 2
Cl[S:2]([CH2:5][CH2:6][CH2:7][NH:8][C:9](=[O:11])[CH3:10])(=[O:4])=[O:3].[OH:19][CH2:20][C:21]([CH3:38])([CH3:37])[C@@H:22]([O:29][Si:30]([CH3:36])([CH3:35])[C:31]([CH3:34])([CH3:33])[CH3:32])/[CH:23]=[CH:24]/[C:25]([O:27][CH3:28])=[O:26]>>[C:9]([NH:8][CH2:7][CH2:6][CH2:5][S:2]([O:19][CH2:20][C:21]([CH3:38])([CH3:37])[C@@H:22]([O:29][Si:30]([CH3:36])([CH3:35])[C:31]([CH3:32])([CH3:34])[CH3:33])/[CH:23]=[CH:24]/[C:25]([O:27][CH3:28])=[O:26])(=[O:4])=[O:3])(=[O:11])[CH3:10], 2
```
The processed data can also be obtained from https://www.dropbox.com/scl/fo/v1rhyes2wvead9dz3x4fb/h?rlkey=nqtst7azldcry3ixnoigmcv3v&dl=0 <br>

## Training
You can run ```train.sh``` files directly to get models of catalyst, solvent, and reagent, we recommend using GPUs for faster training. Corresponding models are also given in Models. <br>
To train model for a particular condition, take solv0 for example, you can run:<br>
```
chemprop_train --data_path ./data/GCN_solv0_data_withN/GCN_data_train.csv --separate_val_path ./data/GCN_solv0_data_withN/GCN_data_val.csv --separate_test_path ./data/GCN_solv0_data_withN/GCN_data_test.csv --dataset_type multiclass --multiclass_num_classes 697 --save_dir ./data/models/GCN_solv0_withN  --reaction --extra_metrics accuracy top3 --epochs 35
```


## Predicting
You can run ```prediction.sh```to get the raw prediction of the test dataset, which will generate a ```raw_prediction``` folder.

# Condition Clusterer
Run ```class_conditions.py``` and it calculates the reaction conditions under each template and clusters them. The files obtained by clustering are also used as candidate libraries for subsequent predictions:
```
python class_conditions.py --Inclusion 0.9 --data_set train --tpl_radius 0
```
```Inclusion```is the tolerance for labels of each category.   It indicates that a label can be selected as a category label only when the number of times it appears is greater than the total number of conditions times Inclusion.<br>
```data_set```data_set is the data set to be collected. In addition to train, test, val, you can also select all.<br>
It will eventually output a json file```conditions_library.json```.

# Using trained model to make predictions
## Files required for the task
To perform the conditional prediction task, we need the following files: <br>
```dMPNN checkpoints```: A trained dMPNN scoring model is used to score each component in a reaction condition.<br>
```conditions_librarys```: The reaction Condition library is extracted from the training data set, and the responses under each template are clustered using Condition Cluster. This file is used to provide a candidate list for conditional predictions.<br>
```Condition keys```: Contains three CSV files for decoding the predicted catalyst, solvent, and reagent.<br>
These files can be found in the ```model``` file.<br>

## Predicting
The file that predicts the input should contain at least two parts, one reflecting the mapped smiles and the other its corresponding template.<br>
```
_id,Mapped Reaction,template
65986967, [CH3:1][O:2][C:3]1=[CH:4][CH:5]=[C:6]([CH:9]=O)[CH:7]=[CH:8]1.[H][C@@:11]([NH2:10])([CH2:12][C:13]1=[CH:14][NH:15][C:16]2=[C:21]1[CH:20]=[CH:19][CH:18]=[CH:17]2)[C:22]([OH:23])=[O:24]>>[CH3:1][O:2][C:3]1=[CH:4][CH:5]=[C:6]([CH:9]2[NH:10][CH:11]([C:22]([OH:23])=[O:24])[CH2:12][C:13]3=[C:14]2[NH:15][C:16]2=[CH:17][CH:18]=[CH:19][CH:20]=[C:21]32)[CH:7]=[CH:8]1, [CH;D3;+0:2]-[NH;D2;+0:3]-[CH;D3;+0:1]-[c;H0;D3;+0:4]>>O=[CH;D2;+0:1].[C@H;D3;+0:2]-[NH2;D1;+0:3].[cH;D2;+0:4]
45394169, [CH3:1][OH:2].[CH3:23][O:22][C:21]1=[CH:20][CH:19]=[C:18]([CH:16]2[NH:17][CH:5]([C:3](O)=[O:4])[CH2:6][C:7]3=[C:8]2[NH:9][C:10]2=[CH:11][CH:12]=[CH:13][CH:14]=[C:15]32)[CH:25]=[CH:24]1>>[CH3:1][O:2][C:3](=[O:4])[CH:5]1[CH2:6][C:7]2=[C:8]([NH:9][C:10]3=[CH:11][CH:12]=[CH:13][CH:14]=[C:15]23)[CH:16]([C:18]2=[CH:19][CH:20]=[C:21]([O:22][CH3:23])[CH:24]=[CH:25]2)[NH:17]1, [C;H0;D3;+0:1]-[O;H0;D2;+0:2]>>O-[C;H0;D3;+0:1].[OH;D1;+0:2]
2729754, [CH3:1][O:2][C:3](=[O:4])[CH:5]1[CH2:6][C:7]2=[C:8]([NH:9][C:10]3=[CH:11][CH:12]=[CH:13][CH:14]=[C:15]23)[CH:16]([C:18]2=[CH:19][CH:20]=[C:21]([O:22][CH3:23])[CH:24]=[CH:25]2)[NH:17]1>>[CH3:1][O:2][C:3](=[O:4])[C:5]1=[CH:6][C:7]2=[C:8]([NH:9][C:10]3=[CH:11][CH:12]=[CH:13][CH:14]=[C:15]23)[C:16]([C:18]2=[CH:19][CH:20]=[C:21]([O:22][CH3:23])[CH:24]=[CH:25]2)=[N:17]1, [c;H0;D3;+0:5]-[c;H0;D3;+0:4]:[n;H0;D2;+0:3]:[c;H0;D3;+0:2]:[cH;D2;+0:1]>>[CH2;D2;+0:1]-[CH;D3;+0:2]-[NH;D2;+0:3]-[CH;D3;+0:4]-[c;H0;D3;+0:5]
28172530, CO[C:22](=[O:23])[C:11]1=[CH:12][C:13]2=[C:14]([NH:15][C:16]3=[CH:17][CH:18]=[CH:19][CH:20]=[C:21]23)[C:9]([C:6]2=[CH:5][CH:4]=[C:3]([O:2][CH3:1])[CH:8]=[CH:7]2)=[N:10]1>>[CH3:1][O:2][C:3]1=[CH:4][CH:5]=[C:6]([C:9]2=[N:10][C:11]([C:22](=[O:23])[NH:24][NH2:25])=[CH:12][C:13]3=[C:14]2[NH:15][C:16]2=[CH:17][CH:18]=[CH:19][CH:20]=[C:21]32)[CH:7]=[CH:8]1, [C;H0;D3;+0:1]>>C-O-[C;H0;D3;+0:1]
```
You can make predictions by running the corresponding ```prediction.py``` directly, the specific use case is as follows:<br>
```
python ./prediction.sh --test_path ./data/test_data.csv --model_path ./data/model/MPNN_models --key_path ./data/model/keys --library_path ./data/model/condition_library --save_path ./data/condition_pred.json
```
You'll end up with a json file that contains predictions for reaction conditions by category:
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



