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
