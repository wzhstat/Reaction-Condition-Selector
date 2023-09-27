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
To get the data for training, you can run ```get_data.py```, this  uses only the CPU and takes about 3 hours to run on a laptop, but before runing the program here are some config you can change:
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




```1976-2016_5+.csv```: It contains the data that the template has more than m reactions.<br>
```classif_by_temp```: The file shows the corresponding reaction for each template.<br>
```all_cat_withoutN.csv```, ```all_solv_withN.csv```, ```all_solv_withoutN.csv```, ```all_solv_withN.csv```, ```all_reag_withoutN.csv```, ```all_reag_withN.csv```: These are statistics on reaction conditions.<br>

```classif_by_temp```and```all_cat_withoutN.csv```, ```all_solv_withN.csv```, ```all_solv_withoutN.csv```, ```all_solv_withN.csv```, ```all_reag_withoutN.csv```, ```all_reag_withN.csv``` can be downloaded directly from data.zip.




