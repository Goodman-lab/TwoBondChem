#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 09:38:29 2023

@author: chingchinglam
"""


import get_descriptors_v2 as gd
from get_label_v4 import label_compo
from scipy.spatial import distance
import numpy as np
import pandas as pd

from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import AllChem

from atidx import AtIdx
from rdkit.Chem.Draw import SimilarityMaps

#################
# get descriptor 

def descriptors_from_mol(mol, descriptor_arr ='2-bond+'):

    react_basic_info=gd.get_basic_info(mol)
    react_ring_info=gd.get_ring_info(mol)

    react_bond=gd.bondtype()
    react_bond.get_bondtype(mol)
    react_bdtype=react_bond.bond_type_ls

    react_connect=gd.connectivity(mol)
    
    #charge_info=gd.get_charge_info(mol)
    #all_info=[react_basic_info,react_ring_info,react_bdtype,react_c,charge_info]

    if descriptor_arr =='2-bond+':
        react_connect.full_connect_info(maxbond=1)
        react_c=react_connect.com_ls

        all_info=[react_basic_info,react_ring_info,react_bdtype,react_c]

    elif descriptor_arr =='2-bond':
        react_connect.full_connect_info(maxbond=1)
        react_c=react_connect.com_ls

        all_info=[react_basic_info,react_c]

    elif descriptor_arr =='1-bond':
        react_connect.full_connect_info(maxbond=0)
        react_c=react_connect.com_ls

        all_info=[react_basic_info,react_c]

    
    atom_len=len(react_basic_info)

    full_descriptor_ls=[]
    for idx in range(0,atom_len):
        sub_ls=[]
        for info in all_info:
            for i in info[idx]:
                sub_ls.append(i)
        full_descriptor_ls.append(sub_ls)
        
    return full_descriptor_ls
        

    
def descriptors_from_rxn(rxn, descriptor_arr ='2-bond+'):


    reaction=AtIdx([rxn])
    reaction.get_info()

    prod_mol=reaction.prod_H
    react_mol=reaction.react_H
    
    
    full_descriptor_react = descriptors_from_mol(react_mol, descriptor_arr = descriptor_arr)
    full_descriptor_prod = descriptors_from_mol(prod_mol, descriptor_arr = descriptor_arr)
    
    
    return full_descriptor_react, full_descriptor_prod



#################
# get label 



def method2(rxns_str, with_sign=False, bl_method = 'rdkit'):

    ## method 2
    ## atom to atom comparison of change in connectivity and strength 
    test = label_compo(rxns_str)
    test.BondChange_info(w = 0.5, method = bl_method)

    label_ls=test.prod_BondChange
    
    label_para=[]
    for ls in label_ls:
        dist=distance.euclidean(0*len(ls), ls)
        
        ## Add a sign to indicate ‘bond breaking or forming’
        ## Negative – Decrease in connectivity
        ## Positive - Increase in connectivity or involved in atom/group transferred process 
        
        
        if with_sign ==True:
            if ls[3]!=0 or (ls[1]!=0 and ls[0]==0):
                dist=dist*-1
        
            if ls[2]!=0:  
                dist=np.absolute(dist)
        
        label_para.append(dist)

    rev_rdict2 = {v: k for k, v in test.rdict2.items()}

    label_v=label_para.copy()

    for i in range(1,len(label_para)+1):
        idx=rev_rdict2.get(i)
        
        label_v[idx] = label_para[i-1]
        
    ## get the correct atom order using dictionary 
    
    prod_idx_ls2=[test.frp_dict.get(idx) for idx in range(0,len(label_v))]
    label_prod=[label_v[idx] for idx in prod_idx_ls2]
    
    ## categorize 
    
    clabel_v = []
    clabel_prod = []
    
    for idx in range(0,len(label_v)):
        
        if label_v[idx] > 0:
            clabel_v.append(1)
        else: 
            clabel_v.append(0)
        
        if label_prod[idx] > 0:
            clabel_prod.append(1)
        else: 
            clabel_prod.append(0)
            
        
    react, prod , frp_dict= get_info_visual(rxns_str)
    
    react_env=get_all_env(react)
    prod_env=get_all_env(prod)

    #clabel_v_v2=label_at_with_same_env(react_env, clabel_v)
    #clabel_prod_v2=label_at_with_same_env(prod_env, clabel_prod)

    clabel_v_v2=label_at_with_same_env(react_env, clabel_v)
    clabel_prod_v2=label_at_with_same_env(prod_env, clabel_prod)

        
    
    return clabel_v_v2 ,clabel_prod_v2


def method3(rxns_str, with_sign=False,bl_method = 'rdkit'):

    ## method 3
    ## atom to atom comparison of change in just the connectivity 
    test = label_compo(rxns_str)
    test.BondChange_info(w = 0.5,method = bl_method)

    label_ls=test.prod_BondChange

    label_para=[]
    

    
    for ls in label_ls:
        
        prod_ls=0*len(ls[2:])
        
        if with_sign ==True:
            ## Add a sign to indicate ‘bond breaking or forming’
            ## Negative – Decrease in connectivity
            ## Positive - Increase in connectivity 
            
            react_ls= [ls[2], -1*ls[3]]
        
        else:
            react_ls= [ls[2], ls[3]]
        
                
        dist=distance.euclidean(prod_ls, react_ls)
        #print(dist)
        
        label_para.append(dist)

    rev_rdict2 = {v: k for k, v in test.rdict2.items()}

    label_v=label_para.copy()

    for i in range(1,len(label_para)+1):
        idx=rev_rdict2.get(i)
        
        label_v[idx] = label_para[i-1]
        
    ## get the correct atom order using dictionary 
    
    prod_idx_ls2=[test.frp_dict.get(idx) for idx in range(0,len(label_v))]
    
    
    label_prod=[label_v[idx] for idx in prod_idx_ls2]
    
    ## categorize 
    
    clabel_v = []
    clabel_prod = []
    
    for idx in range(0,len(label_v)):
        
        if label_v[idx] > 0:
            clabel_v.append(1)
        else: 
            clabel_v.append(0)
        
        if label_prod[idx] > 0:
            clabel_prod.append(1)
        else: 
            clabel_prod.append(0)
            
    react, prod , frp_dict= get_info_visual(rxns_str)
    
    react_env=get_all_env(react)
    prod_env=get_all_env(prod)

    clabel_v_v2=label_at_with_same_env(react_env, clabel_v)
    clabel_prod_v2=label_at_with_same_env(prod_env, clabel_prod)

    
    
    return clabel_v_v2 ,clabel_prod_v2



def multi_path_descriptor_n_label(map_rxn_df, code, method='m3',bl_method='rdkit', descriptor_arr ='2-bond+'):

    sub_df = map_rxn_df[map_rxn_df['code'] == code]

    rxn_ls=sub_df['reaction'].tolist()

    reaction=AtIdx([rxn_ls[0]])
    reaction.get_info()
    
    full_descriptor_react = descriptors_from_mol(reaction.react_H, descriptor_arr = descriptor_arr )

    prod_label_ls=[]
    for r in rxn_ls:
        if method =='m3':
            clabel_v_v2 ,clabel_prod_v2 = method3(r, with_sign=False, bl_method=bl_method)
        elif method == 'm2':
            clabel_v_v2 ,clabel_prod_v2 = method2(r, with_sign=False, bl_method=bl_method)
            
        prod_label_ls.append(clabel_v_v2)


    final_label_ls=[]
    for idx in range(0,len(prod_label_ls[0])):
        v=0
        for ls in prod_label_ls:
            if ls[idx] == 1:
                v=1
        final_label_ls.append(v)

    
    return full_descriptor_react, final_label_ls


    
   
#################
## for visualisation

#from rdkit.Chem.Draw import SimilarityMaps
#SimilarityMaps.GetSimilarityMapFromWeights(prod, label_p, contourLines=10)

def get_info_visual(rxn):
    
    reaction=AtIdx([rxn])
    reaction.get_info()

    prod=reaction.prod_H
    react=reaction.react_H
    
    AllChem.Compute2DCoords(prod)
    AllChem.Compute2DCoords(react)
    
    return react, prod, reaction.frp_dict


###############
## get atom environment via isotopic labelling and conversions to inchi


def list_duplicates(seq):

    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)  
    
    return [[key,locs] for key,locs in tally.items() if len(locs)>0]


def get_at_env(atom_sym, react, atom_ls):

    H_atom_ls=[ls[1] for ls in atom_ls if ls[0] == atom_sym]

    isotope = Chem.Atom(atom_sym)
    isotope.SetIsotope(round(isotope.GetMass())+1)

    du_react_ls=[]
    du_react_inchi_ls=[]

    for atom_index in H_atom_ls:

        new_mol = Chem.Mol(react)


        new_mol.GetAtomWithIdx(atom_index).SetAtomicNum(isotope.GetAtomicNum())
        new_mol.GetAtomWithIdx(atom_index).SetIsotope(isotope.GetIsotope())

        du_react_ls.append(new_mol)

        du_react_inchi_ls.append(Chem.MolToInchi(new_mol) )


    sort_inchi_ls=list_duplicates(du_react_inchi_ls)

    env_ls=[]
    for idx in range(0,len(sort_inchi_ls)):

        new_name=atom_sym+'_'+str(idx)

        idx_ls=sort_inchi_ls[idx][1]

        env_ls.append([new_name,[H_atom_ls[i] for i in idx_ls]])

    
    return env_ls


def get_all_env(react):

    atom_ls=[[atom.GetSymbol(), atom.GetIdx()] for atom in react.GetAtoms()]
    ele_ls= list(set([ls[0] for ls in atom_ls]))

    full_env_ls=[]             

    for ele in ele_ls:
        at_env_ls= get_at_env(ele,react, atom_ls)

        for ls in at_env_ls:
            full_env_ls.append(ls) 

    full_env_idx_ls=[ls[1] for ls in full_env_ls]
    
    return full_env_idx_ls

def label_at_with_same_env(full_env_idx_ls, clabel_v):
    
    positive_idx_ls = [idx for idx in range(0,len(clabel_v)) if clabel_v[idx] > 0]
    same_env_idx_ls = []

    for i in positive_idx_ls:
        for ls in full_env_idx_ls:
            if i in ls:
                same_env_idx_ls.append(ls)

    clabel_v_v2 = clabel_v.copy()

    for ls in same_env_idx_ls:
        max_v=max([clabel_v[i] for i in ls])

        for v in ls:
            clabel_v_v2[v] =max_v
        
    return clabel_v_v2






