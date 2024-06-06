#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 20:38:27 2022

@author: chingchinglam
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import os 
from datetime import date
import sys

def get_atomlist(mol_rdkit_H):

    ## mol = rdkit mol
    ## generate atom list (atom symbol, index, str(symbol+index) and [xyz])
    coor_ls=mol_rdkit_H.GetConformer().GetPositions().tolist()
    atom_list =[]
    for atom,ls in zip(mol_rdkit_H.GetAtoms(), coor_ls):
        atom_list.append([atom.GetIdx(), atom.GetSymbol(), ls, atom.GetAtomicNum()])
    
    atom_info = {atom_list[idx][0]:{'sym':atom_list[idx][1], 'xyz':atom_list[idx][2], 'atno':atom_list[idx][3]} 
                 for idx in range(0,len(atom_list))}
    
    return atom_info


def bond_info_array(mol):
    ## rdkit mol
    ## give a list of bond info incl. atom index and bond order
    bonds_info = [[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble()] for bond in mol.GetBonds()]
    return bonds_info 


def get_bondlength_df(smi_str):
    
    pmol = Chem.MolFromSmiles(smi_str)
    mol = Chem.AddHs(pmol)
    AllChem.EmbedMolecule(mol,randomSeed=0xf00d)
    AllChem.MMFFOptimizeMolecule(mol)


    bond_info=bond_info_array(mol)
    bond_info=bond_info_array(mol)

    atom_info =get_atomlist(mol)

    bond_sym=[]
    bond_length=[]
    for ls in bond_info:
        at1=ls[0]
        at2=ls[1]
        bond=sorted([atom_info.get(at1).get('sym'),atom_info.get(at2).get('sym')])
        bond_sym.append(bond[0]+bond[1])
    
        at1_xyz=np.array(atom_info.get(at1).get('xyz'))
        at2_xyz=np.array(atom_info.get(at2).get('xyz'))
        bond_length.append(np.linalg.norm(at1_xyz - at2_xyz))


    result=pd.DataFrame({'sym':bond_sym, 'length':bond_length})
    result['smi']=smi_str
    
    
    return result


#smi_str_df = pd.read_csv('find_bondlength_smi_12102022.csv')

def DfFromSmi(filename):
    smi_output=[]
    with open(filename) as output:
        for line in output:
            smi_output+= [line]

    smi_ls=[]
    name_ls=[]
    for i in smi_output:
        txt=i.split('\t')
        smi_ls.append(txt[0])
        name_ls.append(txt[1][:-1])
    
    smi_str_df = pd.DataFrame({'name':name_ls,'smi':smi_ls})
    
    return smi_str_df

path = os.getcwd()+'/data/'
today = date.today()
date_str=today.strftime("%d%m%Y")
this_python_script_name= sys.argv[0].split('/')[-1][:-3]


smi_str_df = DfFromSmi('CHEMBL_100k.smi')

smi_str_ls=smi_str_df['smi'].tolist()
bondlength_df_ls=[]

for smi in smi_str_ls:
    try:
        b_df=get_bondlength_df(smi)
        bondlength_df_ls.append(b_df)
    except:
        pass

    
bondlength_df = pd.concat(bondlength_df_ls)

bondlength_df.to_csv(this_python_script_name+'_'+date_str+'.csv')

