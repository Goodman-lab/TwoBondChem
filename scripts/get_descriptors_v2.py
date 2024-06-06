#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 11:03:33 2023

@author: chingchinglam
"""


########### 
## import libraries
## definition works 

from rdkit import Chem
from rdkit.Chem import AllChem

import os

import pandas as pd
import numpy as np
import json

from collections import defaultdict


element_data = {1: {'x': 1, 'y': 1}, 2: {'x': 18, 'y': 1}, 3: {'x': 1, 'y': 2}, 4: {'x': 2, 'y': 2}, 5: {'x': 13, 'y': 2}, 6: {'x': 14, 'y': 2}, 
                7: {'x': 15, 'y': 2}, 8: {'x': 16, 'y': 2}, 9: {'x': 17, 'y': 2}, 10: {'x': 18, 'y': 2}, 11: {'x': 1, 'y': 3}, 12: {'x': 2, 'y': 3}, 
                13: {'x': 13, 'y': 3}, 14: {'x': 14, 'y': 3}, 15: {'x': 15, 'y': 3}, 16: {'x': 16, 'y': 3}, 17: {'x': 17, 'y': 3}, 18: {'x': 18, 'y': 3}, 
                19: {'x': 1, 'y': 4}, 20: {'x': 2, 'y': 4}, 21: {'x': 3, 'y': 4}, 22: {'x': 4, 'y': 4}, 23: {'x': 5, 'y': 4}, 24: {'x': 6, 'y': 4}, 
                25: {'x': 7, 'y': 4}, 26: {'x': 8, 'y': 4}, 27: {'x': 9, 'y': 4}, 28: {'x': 10, 'y': 4}, 29: {'x': 11, 'y': 4}, 30: {'x': 12, 'y': 4}, 
                31: {'x': 13, 'y': 4}, 32: {'x': 14, 'y': 4}, 33: {'x': 15, 'y': 4}, 34: {'x': 16, 'y': 4}, 35: {'x': 17, 'y': 4}, 36: {'x': 18, 'y': 4}, 
                37: {'x': 1, 'y': 5}, 38: {'x': 2, 'y': 5}, 39: {'x': 3, 'y': 5}, 40: {'x': 4, 'y': 5}, 41: {'x': 5, 'y': 5}, 
                42: {'x': 6, 'y': 5}, 43: {'x': 7, 'y': 5}, 44: {'x': 8, 'y': 5}, 45: {'x': 9, 'y': 5}, 46: {'x': 10, 'y': 5}, 47: {'x': 11, 'y': 5}, 
                48: {'x': 12, 'y': 5}, 49: {'x': 13, 'y': 5}, 50: {'x': 14, 'y': 5}, 51: {'x': 15, 'y': 5}, 52: {'x': 16, 'y': 5}, 53: {'x': 17, 'y': 5}, 
                54: {'x': 18, 'y': 5}, 55: {'x': 1, 'y': 6}, 56: {'x': 2, 'y': 6}, 57: {'x': 3, 'y': 9}, 58: {'x': 4, 'y': 9}, 59: {'x': 5, 'y': 9}, 60: {'x': 6, 'y': 9},
                61: {'x': 7, 'y': 9}, 62: {'x': 8, 'y': 9}, 63: {'x': 9, 'y': 9}, 64: {'x': 10, 'y': 9}, 65: {'x': 11, 'y': 9}, 66: {'x': 12, 'y': 9}, 
                67: {'x': 13, 'y': 9}, 68: {'x': 14, 'y': 9}, 69: {'x': 15, 'y': 9}, 70: {'x': 16, 'y': 9}, 71: {'x': 17, 'y': 9}, 72: {'x': 4, 'y': 6}, 
                73: {'x': 5, 'y': 6}, 74: {'x': 6, 'y': 6}, 75: {'x': 7, 'y': 6}, 76: {'x': 8, 'y': 6}, 77: {'x': 9, 'y': 6}, 78: {'x': 10, 'y': 6}, 
                79: {'x': 11, 'y': 6}, 80: {'x': 12, 'y': 6}, 81: {'x': 13, 'y': 6}, 82: {'x': 14, 'y': 6}, 83: {'x': 15, 'y': 6}, 84: {'x': 16, 'y': 6}, 
                85: {'x': 17, 'y': 6}, 86: {'x': 18, 'y': 6}, 87: {'x': 1, 'y': 7}, 88: {'x': 2, 'y': 7}, 89: {'x': 3, 'y': 10}, 90: {'x': 4, 'y': 10}, 
                91: {'x': 5, 'y': 10}, 92: {'x': 6, 'y': 10}, 93: {'x': 7, 'y': 10}, 94: {'x': 8, 'y': 10}, 95: {'x': 9, 'y': 10}, 96: {'x': 10, 'y': 10}, 
                97: {'x': 11, 'y': 10}, 98: {'x': 12, 'y': 10}, 99: {'x': 13, 'y': 10}, 100: {'x': 14, 'y': 10}, 101: {'x': 15, 'y': 10}, 102: {'x': 16, 'y': 10}, 
                103: {'x': 17, 'y': 10}, 104: {'x': 4, 'y': 7}, 105: {'x': 5, 'y': 7}, 106: {'x': 6, 'y': 7}, 107: {'x': 7, 'y': 7}, 108: {'x': 8, 'y': 7}, 
                109: {'x': 9, 'y': 7}, 110: {'x': 10, 'y': 7}, 111: {'x': 11, 'y': 7}, 112: {'x': 12, 'y': 7}, 113: {'x': 13, 'y': 7}, 114: {'x': 14, 'y': 7}, 
                115: {'x': 15, 'y': 7}, 116: {'x': 16, 'y': 7}, 117: {'x': 17, 'y': 7}, 118: {'x': 18, 'y': 7}}






key_element_ls='H, C, N, O, B, F, Cl, Br, Si, P, S'.split(', ')

############# 


def get_basic_info(mol_H, key_element=key_element_ls):
    
    ## input: RDkit mol, with H // Chem.AddHs(mol)
    ## output: [[atomic_number, group, period, isKeyElement],[...],...]
    ## [[14, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],[...],...]
    ## element_data must be in the same script 
    
    atom_list=[]
    for atom in mol_H.GetAtoms():
        atomic_no=atom.GetAtomicNum()
        group = element_data.get(atomic_no).get('x')
        period = element_data.get(atomic_no).get('y')
    
        at_sym = atom.GetSymbol()
        ele_ls=[]
        for j in key_element:
            if at_sym == j:
                ele_ls.append(1)
            else: 
                ele_ls.append(0)
    
        atom_list.append([group, period]+ele_ls)

    return atom_list

##################

def get_ring_info(mol_H):
    
    ## input: RDkit mol, with H // Chem.AddHs(mol)
    ## output: [[IsInRingSize,...],]
    
    
    TF_dict={True:1,False:0}

    ring_size=[3,4,5,6,7]

    ring_info_ls=[]
    for atom in mol_H.GetAtoms():
        ring_info_ls.append([TF_dict.get(atom.IsInRingSize(s)) for s in ring_size]+[TF_dict.get(atom.IsInRing())])

    return ring_info_ls  


#####################
## connectivity

def bond_info_array(mol):
    ## rdkit mol
    ## give a list of bond info incl. atom index and bond order
    bonds_info = [[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble()] for bond in mol.GetBonds()]
    return bonds_info 


def get_BondOrderMatrix(mol_rdkit):
    
    ## get bond order matrix 
    
    am = Chem.GetAdjacencyMatrix(mol_rdkit)
    df = pd.DataFrame(am)
    
    ## float all the element in df
    for col in df.columns: 
        df[col] = df[col].astype('float', errors='ignore')
        
    ## replace the element (1,0) in the df to bond order with bond info list 
    change_list=bond_info_array(mol_rdkit)
    for i in range(0,len(change_list)):
        df[change_list[i][0]][change_list[i][1]]= change_list[i][2]
        df[change_list[i][1]][change_list[i][0]]= change_list[i][2]
    

    return df

def list_duplicates(seq):
    ## from a list of items to [item, [list of indexes at which the item appears in the list]]
    ## e.g. input: [1, 2, 3, 4, 4, 3]
    ## output: [[1, [0]], [2, [1]], [3, [2, 5]], [4, [3, 4]]]

    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)  
    
    return {key:len(locs) for key,locs in tally.items() if len(locs)>0}


def find_neighbour(df_e,atoms_info,key_element,idx):
    
    ## input: bond matrix, atoms_info from get_atomlist(), key_element, atom_idx
    ## output: connectivity_ls: [3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4] neighbour_atom_ls:[1, 12, 13, 14]
    
    at_df=df_e[df_e[idx] != 0]
    ngb_at_ls=list(at_df.index)

    ngb_sym_ls=[atoms_info.get(n).get('sym') for n in ngb_at_ls]
    dup_ngb_sym_dict=list_duplicates(ngb_sym_ls)
    
    connect_by_ele=[]
    for j in key_element:

        v=dup_ngb_sym_dict.get(j)
        if v == None:
            connect_by_ele.append(0)
        else:
            connect_by_ele.append(v)
            
    connect_ls = connect_by_ele+[len(at_df)]
    
    return connect_ls, ngb_at_ls


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



#### write this in class: -- class connectivity:
    
    
class connectivity:
    def __init__(self, mol_H, key_element=key_element_ls):
        
        
        ## optimise the Chem.Mol() at MMFF level 
        self.mol_H=mol_H
        
        AllChem.EmbedMolecule(self.mol_H,randomSeed=0xf00d)
        AllChem.MMFFOptimizeMolecule(self.mol_H)
        
        self.df_e=get_BondOrderMatrix(self.mol_H)
        self.atoms_info = get_atomlist(self.mol_H)
        self.key_element=key_element
        
        
    def get_self_connectivity(self):
        
        ## self connectivity ls
        ## find the neighbours for every atoms in the mol
    
        self.self_connectivity_ls=[]
        self.ngb_ls=[]

        for i in range(0,len(self.df_e)):
            c_ls, n_ls= find_neighbour(self.df_e,self.atoms_info,self.key_element,i)
   
            self.ngb_ls.append(n_ls)
            self.self_connectivity_ls.append(c_ls)
        
        
        ngb_ls2=[]
        for ls in self.ngb_ls:
            sub_ls=[[i,self.atoms_info.get(i).get('atno')] for i in ls]
            sub_ls2 =sorted(sub_ls, key = lambda x: x[1], reverse=True)
            sub_ls3=[at[0] for at in sub_ls2]
    
    
            ngb_ls2.append(sub_ls3)
        
        self.atom_within_xbond = [[idx] for idx in range(0,len(ngb_ls2))]
        
        self.atom_in_xbond = [ngb_ls2[i] for i in range(0,len(ngb_ls2))]
    
    
    def connectivity_above1bonds(self, atom_within_xbond_ls, atom_in_xbond_ls, x_bond):
        
        #print(atom_within_xbond_ls)
        #print(atom_in_xbond_ls)
        
        ngb_ls=[]
        for ls in atom_in_xbond_ls:

    
            sub_ls=[[i,self.atoms_info.get(i).get('atno')] for i in ls]
            sub_ls2 =sorted(sub_ls, key = lambda x: x[1], reverse=True)
            sub_ls3=[at[0] for at in sub_ls2]
            

            ngb_ls.append(sub_ls3)
                
            
        self.atom_within_xbond = [atom_within_xbond_ls[idx]+ngb_ls[idx] for idx in range(0,len(ngb_ls))]
        
        #print(ngb_ls)
        #print(self.atom_within_xbond )
        
        neigh_connect_ls=[] 
        neigh_s_neigh_ls=[]        
        for ls in ngb_ls: 
            sub_ls=[]
            nsub_ls=[]
            for idx in ls:
                c_ls, n_ls= find_neighbour(self.df_e,self.atoms_info,self.key_element,idx)
                
                sub_ls.append([1]+c_ls)
                nsub_ls.extend(n_ls)
        
            neigh_connect_ls.append(sub_ls)
            neigh_s_neigh_ls.append(nsub_ls)
        
        
        self.atom_in_xbond=[]
        for idx in range(0,len(neigh_s_neigh_ls)):
            filt_ls=[]
            for i in neigh_s_neigh_ls[idx]:
                if i not in self.atom_within_xbond[idx]:
                    filt_ls.append(i)
                    
            self.atom_in_xbond.append(filt_ls)
            
        #print(self.atom_in_xbond)
            
        
        neigh_connect_ls2=[]
        for ls in neigh_connect_ls:
            sub_ls=ls
            max_valence=4*(3**(x_bond-1))
            ## assume max valence = 4
            
            miss = max_valence-len(ls)
            if miss<0:
                miss = 0
            for i in range(0,miss):
                ## to ensure that the array has the same size 
                ## append 0 to fill up the spot if the number of neighbour is smaller than 4
                sub_ls.append([0]*(len(self.key_element)+2))
        
            neigh_connect_ls2.append(sub_ls)
            
        
        ## extend all the list together 
        
        neigh_connect_ls3=[]
        for ls in neigh_connect_ls2:
            sub_ls=[]
            for i in ls:
                for j in i:
                    sub_ls.append(j)
            neigh_connect_ls3.append(sub_ls)   
            
            
        self.neigh_connect = neigh_connect_ls3
        
        
        
        
    def full_connect_info(self, maxbond=2):
        
        ## join self connectivity and neighbour connectivity 
        
        self.get_self_connectivity()
        
        self_ls = self.self_connectivity_ls
        
        neigh_ls_ls=[]
        
        
        for i in range(1,maxbond+1):
            self.connectivity_above1bonds(self.atom_within_xbond, self.atom_in_xbond,i)
            neigh_ls_ls.append(self.neigh_connect)
            
        #print(len(neigh_ls_ls))
            
        self.com_ls=[]
        
        for idx in range(0,len(self_ls)):
            
            sub_ls=self_ls[idx]
            
            for ls in neigh_ls_ls:
                sub_ls += ls[idx]
            
            self.com_ls.append(sub_ls)
        
        
        


###############################
## bond strength descriptors


### this function has been introduced to the class bondtype
def get_bondlength_df(mol):
    
    # RDkit mol with H atoms
    
    AllChem.EmbedMolecule(mol)
    
    AllChem.MMFFOptimizeMolecule(mol)


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


    #result=pd.DataFrame({'sym':bond_sym, 'length':bond_length})
    result=[[b,bl] for b,bl in zip(bond_sym,bond_length)]
    
    return result

####

#path = '/'.join(sys.argv[0].split('/')[:-1])+'/'

#bond_class_para = 'bond_classification_01112022.csv'

class bondtype:
    def __init__(self, bond_class_para='bond_classification_01112022.csv'):
        
        ## import and process bond_class_para to get the criteria for classifying the bonds 
        
        path = os.path.dirname(os.path.abspath(__file__))+'/'
        
        BondClass=pd.read_csv(path+bond_class_para)
        bond_ls=BondClass['bond'].tolist()
        boundary_ls=[[0]+json.loads(i)+[3] for i in BondClass['boundary'].tolist()]
        self.label_ls=[]
        for ls,bond in zip(boundary_ls,bond_ls):
            sub_ls=[bond+'_'+str(idx) for idx in range(0,len(ls)-1)]
            self.label_ls.append(sub_ls)


        self.BondClass_dict = {bond:{'boundary':bondary, 'label':label} 
                  for bond,bondary,label in zip(bond_ls,boundary_ls,self.label_ls)}        
        
        
        self.label=[]
        for ls in self.label_ls:
            for i in ls:
                self.label.append(i)
    
    
    def get_bondlength_df(self,mol_H):
        
        ## 1. FF optimisation of the chem.mol
        
        self.mol1=mol_H
        
        AllChem.EmbedMolecule(mol_H,randomSeed=0xf00d)
        #randomSeed=0xf00d
        AllChem.MMFFOptimizeMolecule(mol_H)


        self.bond_info=bond_info_array(mol_H)

        self.atom_info=get_atomlist(mol_H)
        
        ## 2. find the bond length of all the bonds in the chem.mol
        
        bond_sym=[]
        bond_length=[]
        for ls in self.bond_info:
            at1=ls[0]
            at2=ls[1]
            bond=sorted([self.atom_info.get(at1).get('sym'),self.atom_info.get(at2).get('sym')])
            bond_sym.append(bond[0]+bond[1])
        
            at1_xyz=np.array(self.atom_info.get(at1).get('xyz'))
            at2_xyz=np.array(self.atom_info.get(at2).get('xyz'))
            bond_length.append(np.linalg.norm(at1_xyz - at2_xyz))


        #result=pd.DataFrame({'sym':bond_sym, 'length':bond_length})
        self.ff_bond_length_ls=[[b,bl] for b,bl in zip(bond_sym,bond_length)]
        
        
        
    
        
    def get_bondtype(self,mol_H):
        
        
        ## make 
        
        self.get_bondlength_df(mol_H)
        bond_info2=[ls[:2] for ls in self.bond_info]
        
        ## for each atom identify their bond information 
        self.bond_info_toAtom=[]
        for idx in range(0,len(self.atom_info)):
            sub_ls=[]
            for ls in bond_info2:
                if idx in ls:
                    sub_ls.append(bond_info2.index(ls))
            self.bond_info_toAtom.append(sub_ls) 
        
        
        self.classification_ls=[]
        for ls in self.ff_bond_length_ls:
            try:
                #print(ls[0])
                
                ## pd.cut - sort data values into bins
                
                class_data=pd.cut([ls[1]],sorted(self.BondClass_dict.get(ls[0]).get('boundary')), 
                                  labels=self.BondClass_dict.get(ls[0]).get('label'),
                        duplicates='drop').tolist()
                
                #print(class_data)
                
                self.classification_ls.append(class_data[0])
            except:
                self.classification_ls.append('NA_0')
        
        #print(self.classification_ls)
    
        self.bond_info_class=[]
        for ls in self.bond_info_toAtom:
            sub_ls=[self.classification_ls[i] for i in ls]
            self.bond_info_class.append(sub_ls)

        
        ## note down the index that corresponds to the bond type class for which the atom is in 
        
        bond_info_class_idx=[]
        for ls in self.bond_info_class:
            #print(ls)
            sub_ls=[]
            for i in ls:
                #print(i)
                if i != 'NA_0':
                    sub_ls.append(self.label.index(i))
            bond_info_class_idx.append(sub_ls)
    
        #print(bond_info_class_idx)
        
        ## assemble the descriptor list according to the bond_info_class_idx
        
        self.bond_type_ls=[]
        for ls in bond_info_class_idx:
            ini_ls=len(self.label)*[0]
            for i in ls:
                ini_ls[i]+=1
            self.bond_type_ls.append(ini_ls)


#################
## charges 

def get_charge_info(mol_H):
    
    AtFormalCharge_ls=[atom.GetFormalCharge() for atom in mol_H.GetAtoms()]
    
    Chem.rdPartialCharges.ComputeGasteigerCharges(mol_H)
    GasteigerCharge_ls = [atom.GetDoubleProp('_GasteigerCharge') for atom in mol_H.GetAtoms()]
    
    
    overall_charge = Chem.GetFormalCharge(mol_H)
    fragments = Chem.GetMolFrags(mol_H, asMols=True)
    fragments_idx = list(Chem.GetMolFrags(mol_H))
    
    nomol=len(fragments)
    
    frag_charge_ls=[Chem.GetFormalCharge(m) for m in fragments] 
    
    desp_ls=[]
    
    
    for idx in range(0,len(AtFormalCharge_ls)):
    
        at_desp_ls=[nomol, overall_charge]
        at_desp_ls+=[AtFormalCharge_ls[idx]]
        at_desp_ls+=[GasteigerCharge_ls[idx]]
        
        for i in range(0,len(fragments_idx)):
            if idx in fragments_idx[i]:
                self_mol_charge=frag_charge_ls[i]
                other_mol_charge = frag_charge_ls[:i] + frag_charge_ls[i+1:]
                arranged_frag_charge_ls = [self_mol_charge]+other_mol_charge+(3-len(fragments_idx))*[0]
        
        
        at_desp_ls+= arranged_frag_charge_ls
        
    
        desp_ls.append(at_desp_ls)
    
    return desp_ls
    
    




