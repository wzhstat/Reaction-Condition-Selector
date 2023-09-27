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

You can dowload the processed documents from the _data_ file or you can run *get_data.py* program to generate the datasets.<br>

## config settings
Here are some config you can change before generating a dataset:
```
config = {
    'data_path':'./data/grants',
    'save_path':'./data',
    "min_num_covered_rxns_by_rxn_centralized_template":5,
    "min_num_covered_rxns_by_solvent":5,
    "min_num_covered_rxns_by_reagent":5,
    "min_num_covered_rxns_by_catalyst":5,

}
```




