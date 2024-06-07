# TwoBondChem

Research data and scripts supporting 'Every atom counts: Predicting sites of reaction based on chemistry within two bonds'

Directory tree of the TwoBondChem/:

```
TwoBondChem/
├─ bond_strength_descriptor/
│   ├─ bond_pre_classifier_v2.py
│   └─ bondlength_fromsmi_v2.py
├─ data_processing/
│   ├─ competitive_pathway_atom_mapping.py
│   ├─ initial_treatment.py
│   └─ remove_emb_error.py
├─ dataset/
│   ├─ cyclo_data_v2_13012024.csv
│   ├─ da_08012024_vCorr.csv
│   ├─ diels_alder_data_v7_19052024.csv
│   ├─ exam_test.csv
│   ├─ exam_test_28122023.csv
│   ├─ fav_RGD1_13012024_v2.csv
│   ├─ first-year_data_all_27122023.csv
│   ├─ green_04012024.csv
│   ├─ green_ea60_c_23052024.csv
│   └─ rgd_ea60_c_23052024.csv
├─ model_evaluation/
│   ├─ chk_dup_v3.py
│   ├─ random_sample_test.py
│   ├─ random_sample_test_all.py
│   ├─ random_sample_test_Bmix.py
│   ├─ random_sample_test_Cmix.py
│   ├─ save_model.py
│   ├─ save_model2.py
│   ├─ take_one_out_test.py
│   └─ train_size_test.py
├─ model_test_2.ipynb
├─ React.py
├─ README.md
├─ RF_model_yr1_28122023.sav
└─ scripts/
    ├─ atidx.py
    ├─ bond_classification_01112022.csv
    ├─ get_descriptors_v2.py
    ├─ get_label_v4.py
    ├─ model_training_v3.py
    └─ training_prepare_v5.py

*Archive files have not been included 

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

We filtered and processed the RGD1 and the Green datasets to ensure that the reactions are thermodynamically favourable (ie ΔHr < 0 kcal mol-1) with a low kinetic barrier (ie activation energy, EA < 40 kcal mol-1). This gives the ‘fav_RGD1_13012024_v2.csv’ and ‘green_04012024.csv’. To study the effect of varying the EA cut-off, we also filtered and processed the RGD1 and the Green datasets with an EA cut-off of 60 kcal mol and ΔHr < 0 kcal mol-1 – ‘rgd_ea60_c_23052024.csv’ and ‘green_ea60_c_23052024.csv’.

data_processing/competitive_pathway_atom_mapping.py input and output example files: ‘exam_test.csv’ and ‘exam_test_28122023.csv’


## C. Scripts 


React_v2.py: code for using the random forest model

model_test_3.ipynb: Jupyter note on how to use the ‘React.py’ script 


Saved models: 



### scripts/ folder: 

atidx.py: contains the functions for conducting atom-to-atom mapping and formatting the reaction SMILES strings for competitive pathways 


get_descriptors_v2.py: contains functions for generating the atomistic descriptor components 

-   Associated file: ‘bond_classification_01112022.csv’ – the parameters for generating the bond strength descriptors 


get_label_v4.py: contains functions for generating the atomistic label 

training_prepare_v5.py: compiles functions in ‘get_descriptors_v2.py’ and ‘get_label_v4.py’ to generate the descriptor arrays and labels for atoms in a set of reactants

model_training_v3.py: functions for training and evaluating the model 


### data_processing/ folder: 


initial_treatment.py: group reactions with the same reactants together vis InChI strings; for the Diels-Alder dataset only, reactions with placeholder atoms are also removed. Reactions with hypervalent molecules lead to errors in the model training and testing. They are removed manually from the Diels-Alder dataset.


remove_emb_error.py: filter reactions with 3D structure errors


competitive_pathway_atom_mapping.py: This script executes the functions in atidx.py to conduct atom-to-atom mapping of individual reactions and between the competitive pathways. This process was conducted only on the first-year and Diels-Alder datasets. 

-   Associated files in the dataset/ folder: ‘exam_test.csv’ and ‘exam_test_28122023.csv’ – the input and output csv file from executing the script are provided for illustrations 


### bond_strength_descriptor/ folder: 


bondlength_fromsmi_v2.py: gather bond length data from SMILES strings 

bond_pre_classifier_v2.py: using the bond length data to generate bond strength descriptor classification criteria, ie the ‘bond_classification_01112022.csv’ file


### model_evaluation/ folder: 


chk_dup_v3.py: check if there is an overlap between datasets


random_sample_test.py, random_sample_test_all.py, random_sample_test_Bmix.py, random_sample_test_Cmix.py: scripts for carrying out random sampling test 


take_one_out_test.py: the script for carrying out the take-one-out cross-validation test 


test_size_test.py: The script for investigating the effect of varying the testing or training dataset size 


save_model.py: executing this script to train models using reactions from a single dataset and return the trained models in .sav files

save_model2.py: executing this script to train models using reactions from multiple datasets and return the trained models in .sav files




