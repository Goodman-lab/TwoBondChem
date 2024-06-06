#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 17 13:30:18 2024

@author: chingchinglam
"""

from datetime import date
import sys
import os
sys.path.insert(1, os.getcwd()+'/scripts/')

import pandas as pd 
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem

def list_duplicates(seq):
    ## from a list of items to [item, [list of indexes at which the item appears in the list]]
    ## e.g. input: [1, 2, 3, 4, 4, 3]
    ## output: [[1, [0]], [2, [1]], [3, [2, 5]], [4, [3, 4]]]

    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)  
    
    return [[locs,key] for key,locs in tally.items() if len(locs)>0]

def parah_code(given_list):

    first_value=given_list[0]
    first_code=0

    new_code_ls=[0]

    for v in given_list[1:]:
        
        if v != first_value:
            first_value=v
            first_code += 1

            new_code_ls.append(first_code)

        else:
            new_code_ls.append(first_code)
    
    return new_code_ls

full_data= pd.read_csv('rgd_ea60.csv')
rxn_ls=full_data['reaction'].tolist()

react_ls=[]
prod_ls=[]
for rxn in rxn_ls:
    sub_ls=rxn.split('>>')
    react_ls.append(sub_ls[0])
    prod_ls.append(sub_ls[1])

ps = Chem.SmilesParserParams()
ps.removeHs = False

error_idx_ls=[]

for idx in range(0,len(react_ls)):
    try:
        mol=Chem.MolFromSmiles(react_ls[idx],ps)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        AllChem.EmbedMolecule(mol,randomSeed=0xf00d)
        AllChem.MMFFOptimizeMolecule(mol)
        coor_ls=mol.GetConformer().GetPositions().tolist()
    except:
        error_idx_ls.append(idx)

for idx in range(0,len(prod_ls)):
    try:
        mol=Chem.MolFromSmiles(prod_ls[idx],ps)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        AllChem.EmbedMolecule(mol,randomSeed=0xf00d)
        AllChem.MMFFOptimizeMolecule(mol)
        coor_ls=mol.GetConformer().GetPositions().tolist()
    except:
        error_idx_ls.append(idx)

print(error_idx_ls)        

        
rxn_idx_ls=[idx for idx in range(len(rxn_ls))]
ef_rxn_idx_ls=[i for i in rxn_idx_ls if i not in error_idx_ls]

sample_df=full_data.loc[full_data['idx'].isin(ef_rxn_idx_ls)].copy()

sample_df['code'] = parah_code(sample_df['code'].tolist())
sample_df['idx'] = [i for i in range(0,len(sample_df))]



file_name = 'rgd_ea60_c'
today = date.today()
date_str=today.strftime("%d%m%Y")
sample_df.to_csv(file_name+'_'+date_str+'.csv' )
