#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 13:46:58 2023

@author: chingchinglam


v2 -- update class AtIdx
v3 -- remove class AtIdx; class AtIdx is now in atidx.py
v4 -- update on the label for change in bond strength 

"""

import numpy as np
from rdkit import Chem
#from rxnmapper import RXNMapper
import get_descriptors_v2 as gd
from atidx import AtIdx
    
            

##############################
###########


def get_bond_as_idx(rdkit_H,bond_info_class,bond_info_toAtom, method = 'bond_length'):
    ## get bond info of a rdkit mol in [[1, 2, 'CC', '1']...] form 
    
    react_bonds=[sorted(i[:2]) for i in gd.bond_info_array(rdkit_H)]
    react_H_at = [at for at in rdkit_H.GetAtoms()]

    react_bonds1=[]
    for ls in react_bonds:
        react_bonds1.append([react_H_at[i].GetAtomMapNum() for i in ls])
        
    react_bonds2 = [sorted(i) for i in react_bonds1]
    
    react_b=bond_info_class

    bond_as_atidx=[]
    
    
    if method == 'bond_length':
    
        for ls,bls in zip(bond_info_toAtom,react_b):
        
            for idx,b in zip(ls,bls):
                bond=b 
            
                if bond == 'NA_0':
                    ## label as NA_0 if it doesn;t appear in the bondtype.BondClass_dict
                
                    sort_bond=react_bonds[idx]
                    sort_bond2=sorted([react_H_at[sort_bond[0]].GetSymbol(),react_H_at[sort_bond[1]].GetSymbol()])
                    bond = ''.join(sort_bond2)+'_0'
            
            
                bond_as_atidx.append(react_bonds2[idx]+bond.split('_'))
                
    elif method == 'rdkit':
        
        for ls,bls in zip(bond_info_toAtom,react_b):
        
            for idx,b in zip(ls,bls):
                
                sort_bond=react_bonds[idx]
                sort_bond2=sorted([react_H_at[sort_bond[0]].GetSymbol(),react_H_at[sort_bond[1]].GetSymbol()])
                
                rdkitbond=rdkit_H.GetBondBetweenAtoms(react_H_at[sort_bond[0]].GetIdx(), react_H_at[sort_bond[1]].GetIdx())
                bondtype=rdkitbond.GetBondTypeAsDouble()
                
                bond = ''.join(sort_bond2)+'_'+str(bondtype)
                
                bond_as_atidx.append(react_bonds2[idx]+bond.split('_'))
                
                
    
    
    bond_as_atidx2 = []

    for i in bond_as_atidx:
        if i not in bond_as_atidx2:
            bond_as_atidx2.append(i)


    bond_as_atidx3 = []
    for idx in range(1,len(react_H_at)+1):
        sub_ls=[]
        for b in bond_as_atidx2:
            if idx in b: 
                sub_ls.append(b)
        bond_as_atidx3.append(sub_ls)

    
    
    return bond_as_atidx3



         
class label_compo:
    def __init__(self,rxn):
        
        ## perform atom to atom mapping 
        
        test=AtIdx([rxn])
        test.get_info()
        
        self.react_H = test.react_H
        self.prod_H = test.prod_H
        
        self.fpr_dict = test.fpr_dict
        self.rdict2 = test.rdict2
        self.frp_dict = test.frp_dict
        
        
        #self.confidence = test.confidence
    
    def connectivity_info(self):
        
        ## using class connectivity from get_descriptor to get the neighbour connectivity info
        
        connect_info=gd.connectivity(self.react_H)
        connect_info.full_connect_info()
        connect_info2=gd.connectivity(self.prod_H)
        connect_info2.full_connect_info()
        
        self.rea_connect_info=connect_info.self_connectivity_ls
        prod_connect_info1=connect_info2.self_connectivity_ls
        
        self.prod_connect_info=[]
        for i in range(0,len(prod_connect_info1)):
            self.prod_connect_info.append(prod_connect_info1[self.fpr_dict.get(i)])
            
            
    def strength_info(self):
        
        ## using class bondtype from get_descriptor to get the bond strength info
        
        bond_info=gd.bondtype()
        bond_info.get_bondtype(self.react_H)

        bond_info2=gd.bondtype()
        bond_info2.get_bondtype(self.prod_H)
        
        self.rea_bond_info=bond_info.bond_type_ls
        prod_bond_info1=bond_info2.bond_type_ls
        
        self.prod_bond_info=[]
        for i in range(0,len(prod_bond_info1)):
            self.prod_bond_info.append(prod_bond_info1[self.fpr_dict.get(i)])
        
        
    def BondChange_info(self, w = 0.5, method = 'rdkit'):
        
        
        ## atom to atom comparison of change in connectivity 
        ## bond strength classification of the bonds that are connected to an atom  

        
        bond_info=gd.bondtype()
        bond_info.get_bondtype(self.react_H)

        bond_info2=gd.bondtype()
        bond_info2.get_bondtype(self.prod_H)
        
        
        ## use get_bond_as_idx to get bond info in [[1, 2, 'CC', '1']...] form
        
        
        react_bond_as_atidx = get_bond_as_idx(self.react_H,bond_info.bond_info_class,bond_info.bond_info_toAtom,method=method)
        prod_bond_as_atidx = get_bond_as_idx(self.prod_H,bond_info2.bond_info_class,bond_info2.bond_info_toAtom,method=method)
        
        
        ## use the bond_as_atidx to categorise the change:
        ## Increase in bond strength: A
        ## Decrease in bond strength: B

        ## Increase in connectivity: C
        ## Decrease in connectivity: D

        
        bond_strength_only=[]
        self.bond_change_react=[]
        self.bond_change_prod=[]
        for i_ls,j_ls in zip(react_bond_as_atidx,prod_bond_as_atidx):
            
            #print(i_ls)
            #print(j_ls)
            #print()

            sub_ls1=[i for i in i_ls]
            #print(sub_ls1)

            sub_ls2=[j for j in j_ls ]
    
            sub_bond_strength_only=[]
            sub_bond_strength_only2=[]
            sub_bond_strength_only3=[]
            for b1 in sub_ls1:
                for b2 in sub_ls2: 
                    if b1[:3] == b2[:3]:
                        sub_bond_strength_only.append(b1[:3]+[float(b1[-1])]+[float(b2[-1])])
                        sub_bond_strength_only2.append(b1)
                        sub_bond_strength_only3.append(b2)


    
            bond_strength_only.append(sub_bond_strength_only)
    
            sub_bond_change_react=[i for i in sub_ls1 if i not in sub_bond_strength_only2]
            sub_bond_change_prod=[i for i in sub_ls2 if i not in sub_bond_strength_only3]
    
            self.bond_change_react.append(sub_bond_change_react)
            self.bond_change_prod.append(sub_bond_change_prod)
            
        #print(bond_strength_only)
        self.bond_strength_up=[]
        self.bond_strength_down=[]
        for ls in bond_strength_only:
            
            
            sub_ls=[]
            for i in ls:
                dl=np.absolute(i[3]-i[4])
                
                if method == 'bond_length':
                    if dl>1 and i[2]!='CH':
                    ## change in CH bond strength does not count - less noisy 
                        sub_ls.append(i)
                else:
                    if dl != 0:
                        sub_ls.append(i)
                        

            
            up_ls=[]
            down_ls=[]
            for j in sub_ls:
                if j[3]>j[4]:
                    up_ls.append(j)
                else:
                    down_ls.append(j)
                
            
        
            self.bond_strength_up.append(up_ls)
            self.bond_strength_down.append(down_ls)
            
            
            
        
        len_bond_strength_up=[len(ls) for ls in self.bond_strength_up]
        len_bond_strength_down=[len(ls) for ls in self.bond_strength_down]
        len_bond_change_react=[len(ls) for ls in self.bond_change_react]
        len_bond_change_prod=[len(ls) for ls in self.bond_change_prod]
        
        
        ## summarising the change in bond strength and connectivity 
        
        self.prod_BondChange=[]
        for idx in range(0,len(len_bond_strength_up)):
            a=len_bond_strength_up[idx]
            b=len_bond_strength_down[idx]
            if a>0:
                a=w
            if b>0:
                b=w
            sub_ls=[a,b,len_bond_change_prod[idx],len_bond_change_react[idx]]
    
            self.prod_BondChange.append(sub_ls)
            
         
        self.react_BondChange=[[0,0,0,0]*len(self.prod_BondChange)]




