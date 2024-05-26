# Reaction-Condition-Selector
This repository contains a reaction condition selector.

# Python Version
python 3.10

# Install Requirements
To run the reaction condition selector you need to make sure you have installed anaconda. The version about ```pytorch``` and ```cudatoolkit``` should be depended on your machine.
```
conda create -n Cluster_Predictor python=3.10 \
conda activate Cluster_Predictor \
pip3 install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cu116 \
pip install pandas==2.1.0 \
pip install chemprop==1.6.1 \
conda install rdkit=2023.3.3 -c rdkit \
```
We also provide the corresponding environment file, You can reproduce the environment directly through the provided .yaml file.<br>
```
conda env create -f environment.yaml
```

# Training D-MPNN Model
## Step 0
Please first download the Cluster Predictor program from GitHub. It contains some of the data,the pre-trained model and the code needed to replicate the work. : <br>
```
git clone https://github.com/wzhstat/Reaction-Condition-Selector.git  \
cd ./Reaction-Condition-Selector 
```
## Step 1 Data Pre-processing
Before training, make sure you have the data file, which should contain two parts: the csv files```data_train.csv```, ```data_test.csv```, ``` data_val.csv``` and the ```labels``` folder. The labels file can be downloaded directly from ```data```, the csv files can be download in https://www.dropbox.com/scl/fo/v1rhyes2wvead9dz3x4fb/h?rlkey=nqtst7azldcry3ixnoigmcv3v&dl=0. You should put them like this: <br>
```
./data/data_train.csv
./data/data_test.csv
./data/data_val.csv
./data/labels
```
You can run the ```preprocessing.py``` file to get the preprocessed data. <br>
```
python preprocessing.py
```


## Step 2 Training
You can run ```train.sh``` files directly to get models of catalyst, solvents, and reagents, we recommend using GPUs for faster training. <br>
To train model for a particular condition constitute, take solv0 for example, you can run:<br>
```
chemprop_train --target_columns solv0 --data_path ./data/MPNN_data/GCN_data_train.csv  --separate_val_path ./data/MPNN_data/GCN_data_val.csv --separate_test_path ./data/MPNN_data/GCN_data_test.csv --dataset_type multiclass --multiclass_num_classes 542 --save_dir ./models/GCN_solv0  --reaction --extra_metrics accuracy --epochs 35
```
You end up with a models folder that contains the trained D-MPNN model. You can also find trained models directly from ```models```. <br>


# Using Trained Model to Make Predictions

## Step 1 Construction of Reaction Condition Library
Before making a prediction, run ```BuildConditionLibrary.sh``` to build the three types of reaction condition libraries needed for prediction, namely the r1 library, the r0 library, and the r0* library. 
```
chmod u+x ./BuildConditionLibrary.sh
./BuildConditionLibrary.sh
```
This will generate a ```condition_library``` folder under the data folder, which contains the following libraries:<br>
```
./data/condition_library/condition_library_r1.json.gz
./data/condition_library/condition_library_r0.json.gz
./data/condition_library/condition_library_r0_1.json.gz
./data/condition_library/classed_conditions_library_r1.json.gz
./data/condition_library/classed_conditions_library_r0.json.gz
./data/condition_library/classed_conditions_library_r0_1.json.gz
```
You can also run ```get_condition_library.py``` directly to get a specific library of reaction conditions.
```
python get_condition_library.py --Inclusion 0.8 --data_set train --tpl_radius 0
```
```Inclusion```is the tolerance for labels of each category.   It indicates that a label can be selected as a category label only when the number of times it appears is greater than the total number of conditions times Inclusion.<br>
```data_set``` is the data set to be collected. <br>
```--tpl_radius``` is the radius of templates used for categorization.<br>

## Step 2 Prediction
You can run ```make_predictions.sh``` to complete the model predictions, which generates a prediction folder containing predictions with and without clustering.<br>
```
chmod u+x ./make_predictions.sh
./make_predictions.sh
```
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
This may take a while, and the calculated result will be printed out：
```
Top-1 Accuracy: 0.44624893250539177
Top-3 Accuracy: 0.6360386179744381
Top-5 Accuracy: 0.706326805332407
Top-10 Accuracy: 0.7854645018599737
Cluster Top-1 Accuracy: 0.6553620796966144
Cluster Top-3 Accuracy: 0.8564770608073878
Cluster Top-5 Accuracy: 0.9074066032683428
Cluster Top-10 Accuracy: 0.9645355855660255
```

# Predictions of conditions in actual drug synthesis routes
## Step 1 Obtaining route data
Our actual route data was artificially extracted from the article, and the processed ```JMC_data.csv``` file can be found in the ```data``` folder.<Br>

## Step 2 Prediction
You can run ```make_JMC_predictions.sh``` to complete the model predictions, which generates a prediction folder containing predictions with and without clustering.<br>
```
chmod u+x ./make_JMC_predictions.sh
./make_JMC_predictions.sh
```

## Setp 3 Calculate Accuracy
Please run Score.py to calculate the accuracy of the model predictions.
```
python Score.py --data_path ./data/JMC_data.csv --pred_path ./data/JMC_prediction
```




