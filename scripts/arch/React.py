#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 10:09:27 2023

@author: chingchinglam
"""

import os 
import sys 
sys.path.insert(1, os.getcwd()+'/scripts/')

import pickle
import training_prepare_v4 as tp 
from rdkit import Chem
from rdkit.Chem import AllChem


class yr1_reactivity:
    
    def __init__(self, model_file='RF_model_yr1_28122023.sav'):
        
        self.model = pickle.load(open(model_file, 'rb'))
        
    
    def make_prediction(self,input_mol):
        
        mol=Chem.MolFromSmiles(input_mol)
        mol_H=Chem.AddHs(mol)
        
        #generate the descriptor
        descriptors_ls=tp.descriptors_from_mol(mol_H)
        #make the prediction
        self.pred=self.model.predict(descriptors_ls)
        
        #for visualisation
        self.mol_H_visual=Chem.Mol(mol_H)
        AllChem.Compute2DCoords(self.mol_H_visual)