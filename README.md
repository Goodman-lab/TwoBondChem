# TwoBondChem

Research data and scripts supporting 'Every atom counts: Predicting sites of reaction based on chemistry within two bonds'

Directory tree of the TwoBondChem/:

```
TwoBondChem/
├─ dataset/
│   ├─ cyclo_data_v2_13012024.csv
│   ├─ da_08012024_vCorr.csv
│   ├─ da_17012024.csv
│   ├─ exam_test_28122023.csv
│   ├─ fav_RGD1_13012024_v2.csv
│   ├─ first-year_data_all_27122023.csv
│   └─ green_04012024.csv
├─ model_test_2.ipynb
├─ React.py
├─ RF_model_yr1_28122023.sav
└─ scripts/
    ├─ atidx.py
    ├─ bond_classificaation_01112022.csv
    ├─ competitive_pathway_atom_mapping.py
    ├─ exam_test.csv
    ├─ exam_test_28122023.csv
    ├─ get_descriptors_v2.py
    ├─ get_label_v4.py
    ├─ model_training_v2.py
    ├─ save_model.py
    └─ training_prepare_v4.py

```
## A. Environment 

The scripts in this project are written in Python under the following environment:

-	The Python (3.8.12) Standard Library
-	Pandas (1.1.5)
-	Numpy (1.19.2) 
-	Sklearn (0.24.2/1.0.1)
-	RDkit (2021.9.4)
-	Scipy (1.4.1)
-	Rxnmapper (0.1.4)

## B. Data 

dataset/ folder:

All the datasets (Type A: the first-year reactions + three reactions from Part 1A exam at Cambridge; Type B: [3+2] cycloaddition1 and Diels-Alder reaction dataset; Type C: the Reaction Graph Depth 1 (RGD1)) and the Green dataset). are available in the dataset/ folder as csv files. The reaction data is processed and formatted in the same style for this investigation. The atom-to-atom mapping numbering in the reaction SMILES is consistent for competitive pathways with the same reactants. There are three columns in each csv file: 

-	idx: index of the reaction
-	code: reactions with the same code are competitive pathways (ie having the same reactants)
-	reaction: mapped reaction SMILES


The RGD1 dataset csv file also has an ‘Rind’ column, which corresponds to the index in the ‘reaction’ column in the original datafile: 
https://figshare.com/articles/dataset/model_reaction_database/21066901?file=40272727). 


‘da_08012024_vCorr.csv’ contains 100 Diels-Alder reactions, where the atom-to-atom mapping errors have been picked out and corrected manually.5 Before the correction, reactions with index = 2, 30, 31, 49, 66, 70, 88 and 99 contain errors. 


Three out of 147 reactions in the first-year reaction dataset (‘first-year_data_all_27122023.csv’) had mapping errors, which were subsequently corrected manually. Before the correction, reactions with index = 50, 74 and 124 contain errors. Reactions with a code greater than 77 belong to the testing data set. 


## C. Scripts 


RF_model_yr1_28122023.sav: A saved random forest model trained from all the first-year reactions

React.py: code for using the random forest model

model_test_2.ipynb: Jupyter note on how to use the ‘React.py’ script 


### The files in the scripts/ folder: 

atidx.py: contains the functions for conducting atom-to-atom mapping and formatting the reaction SMILES strings for competitive pathways 
competitive_pathway_atom_mapping.py: this script executes the functions in atidx.py 

-	Associated files: ‘exam_test.csv’ and ‘exam_test_28122023.csv’ – the input and output csv file from executing the script are provided for illustrations 



get_descriptors_v2.py: contains functions for generating the atomistic descriptor components 

-	Associated file: ‘bond_classificaation_01112022.csv’ – the parameters for generating the bond strength descriptors 


    
get_label_v4.py: contains functions for generating the atomistic label 

training_prepare_v4.py: compiles functions in ‘get_descriptors_v2.py’ and ‘get_label_v4.py’ to generate the descriptor arrays and labels for atoms in a set of reactants 

model_training_v2.py: functions for training and evaluating the model 
save_model.py: executing this script generates the ‘RF_model_yr1_28122023.sav’ file given the ‘first-year_data_all_27122023.csv’ files. 









