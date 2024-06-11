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
import training_prepare_v5 as tp 
from rdkit import Chem
from rdkit.Chem import AllChem


class reactivity:
    
    def __init__(self, data='yr1', model_file='RF_model_yr1_28122023.sav'):

        if data == 'yr1':
            model_file= 'RF_model_yr1_28122023.sav'
        elif data == 'da':
            model_file= 'RF_model_da_11062024.sav'
        elif data == 'cyclo':
            model_file = 'RF_model_cyclo_11062024.sav'
        elif data == 'green':
            model_file = 'RF_model_green_11062024.sav'
        elif data == 'rgd':
            model_file = 'RF_model_rgd_11062024.sav'
        elif data == 'Bmix':
            model_file = 'RF_model_Bmix_11062024.sav'
        elif data == 'Cmix':
            model_file = 'RF_model_Cmix_11062024.sav'
        elif data == 'all':
            model_file = 'RF_model_all_11062024.sav'



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