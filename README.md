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
>This file is used to process the reaction condition data of USPTO1976-2016, and requires the 1976_Sep2016_USPTOgrants_cm file, which can be downloaded from the following address: https://figshare.com/articles/dataset/Chemical_reactions_from_US_patents_1976-Sep2016_/5104873 <br>

## Config Settings
Here are some config you can change before generating a dataset:
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
* ```min_num_covered_rxns_by_catalyst```: This is the minimum number of times a certain catalyst appears

## Output Files
The following files are generated when you run the program: <br>
```1976-2016.csv```: It contains information extracted directly from the USPTO data set, includes *id*, *reaction*, *products*, *reactants*, *reagent*, *catalyst*, *solvent*. Reagent, catalyst and solvent are stored in list form.<br>
```1976-2016_5+.csv```: It contains the data that the template has more than m reactions.<br>
```classif_by_temp```: The file shows the corresponding reaction for each template.<br>




