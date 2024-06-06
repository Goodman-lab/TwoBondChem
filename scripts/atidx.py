#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 12:57:11 2023

@author: chingchinglam

based on class AtIdx from get_label_v2

"""


from rdkit import Chem
import pandas as pd 
from rxnmapper import RXNMapper


def flatten(iterable):
    out = []
    for i in iterable:
        if hasattr(i,'__iter__'):
            out.extend(flatten(i))
        else:
            out.append(i)
    return out


class AtIdx:
    def __init__(self,rxn):
        
        
        ## test if atom to atom mapping is required
        
        self.use_RXNMapper = False
        
        test_mol=Chem.MolFromSmiles(rxn[0].split('>>')[0])
        AtomMapNum_test_mol=[at.GetAtomMapNum() for at in test_mol.GetAtoms()]
        
        if 0 in  AtomMapNum_test_mol:
            self.use_RXNMapper = True
        
        if self.use_RXNMapper == True:
            
            ## atom to atom mapping with RXNMapper
            rxn_mapper = RXNMapper()
            results = rxn_mapper.get_attention_guided_atom_maps(rxn)
            result1=results[0].get('mapped_rxn').split('>>')
            
            self.confidence = results[0].get('confidence')
            
            ## import the result as rdkit mol 
            mol_ls =[ Chem.MolFromSmiles(i) for i in result1]
            
        else:
            
            ps = Chem.SmilesParserParams()
            ps.removeHs = False
            
            map_smi=rxn[0].split('>>')
            
            ## import as rdkit mol 
            mol_ls =[Chem.MolFromSmiles(i,ps) for i in map_smi]
            
        
        
        self.react=mol_ls[0]
        self.prod=mol_ls[1]
        
        Chem.Kekulize(self.react, clearAromaticFlags=True)
        Chem.Kekulize(self.prod, clearAromaticFlags=True)
        


    def get_Hidx(self,mol):
        
        ## process the raw rdkit mol 
        ## obtain indexing info and add H atoms 
        
        len_heavy_at=len(mol.GetAtoms())
        self.at_mapidx=[at.GetAtomMapNum() for at in mol.GetAtoms()]
        
        
        self.at_dict={ i:{'Midx':j} for i,j in zip(range(0,len_heavy_at),self.at_mapidx)}
        
        self.mdict={i:j for i,j in zip(range(0,len(self.at_mapidx)), self.at_mapidx)}
        

        self.mol_H = Chem.AddHs(mol)
    
        mol_H_at = [at for at in self.mol_H.GetAtoms()]
    
        for at in mol_H_at[len_heavy_at:]:   
            at.SetAtomMapNum(at.GetIdx()+1)
    

        H_indx_ls=[]
        for atom in self.mol_H.GetAtoms():   
        
            ngh_ls=atom.GetNeighbors()
            H_indx=[]
            for n in ngh_ls:
                if n.GetSymbol() == 'H':
                    H_indx.append(n.GetIdx()+1)
            
    
            H_indx_ls.append(H_indx)

    
        self.H_indx_ls2=H_indx_ls[:len_heavy_at]
        
        
        for i,j in zip(range(0,len_heavy_at),self.H_indx_ls2):
            self.at_dict[i]['H_ls']=j
            
        
    def process_rxn(self):
        
        ## atom to atom mapping for H atoms 
        ## creating dictionary for mataching indeces 
        ## only applicable if there is only 1 H involved in changes in the connectivity
        
        ## self.get_Hidx() on the reactant and product system 
        self.get_Hidx(self.react)
        
        #print(self.at_mapidx)
        
        r_idx=self.at_mapidx
        rdict=self.mdict
        
        rH_idx=self.H_indx_ls2
        
        self.react_H=self.mol_H
        
        self.get_Hidx(self.prod)
        
        p_idx=self.at_mapidx
        pdict=self.mdict
        pH_idx=self.H_indx_ls2
        self.prod_H=self.mol_H

        #print(rdict)
        #print(pdict)

        ## rdkit idx reactant : product dictionary; heavy atoms only
        pr_dict={}
        rp_dict={}
        for i in range(0,len(r_idx)):
            
            for j in range(0,len(r_idx)):
                #ßprint(i, j)
                if rdict.get(i)==pdict.get(j):
                    #print('x ',i,j)
                    pr_dict[i]=j
                    rp_dict[j]=i
                    
        
        #print(rp_dict)
        #print(len(r_idx))

        
        ## using the rdkit idx reactant : product dictionary to sort the rH_idx (the reactant H atom list)
        sort_rH_idx=[rH_idx[rp_dict.get(i)] for i in range(0,len(r_idx))]
        
        
        
        ###########
        ## the first filter -- 1. chk if there's a H transfer; 
        ## 2. correct the idx of the transfer (i.e. ensure that the idx is the same in the reactant and product system)
        

        match_ls=[]
        for i,j in zip(sort_rH_idx,pH_idx):
            if i != j and len(i) != len(j):
                match_ls.append([False,False])
            elif i != j and len(i) == len(j):
                match_ls.append([False,True])
            else:
                match_ls.append([True,True])
        
        #print(sort_rH_idx)
        #print()
        #print(pH_idx)
        #print(match_ls)
        
        shift_at=[i for i in range(0,len(match_ls)) if match_ls[i] == [False,False]]
        aff_at=[]
        if len(shift_at)!=0:
            aff_at=[i for i in range(0,len(match_ls)) if match_ls[i] == [False,True]]
        
        #aff_at=[i for i in range(aff_at1[0],aff_at1[-1])]
        naff_at=[i for i in range(0,len(match_ls)) if match_ls[i] == [True,True]]

        naff_at_ls=[pH_idx[i] for i in naff_at]
        naff_at_ls2=[]
        for ls in naff_at_ls:
            for i in ls:
                naff_at_ls2.append(i)
                
        #print(shift_at)
        #print(aff_at)
        #print(naff_at)
        
        ## identify the transfer type
        ## L2R = from reactant to product: the proton has moved from the front of the sort_rH_idx list to the back
        ## R2L = from reactant to product: the proton has moved from the back of the sort_rH_idx list to the front 
        ## print sort_rH_idx and pH_idx to check 
        
        if len(shift_at) !=0:
                
            transfer_type='L2R'
            tt_v=-1
        
            if len(pH_idx[shift_at[0]])>len(sort_rH_idx[shift_at[0]]):
                
                #print(len(sort_rH_idx[shift_at[0]]))
                #print(len(pH_idx[shift_at[0]]))
            #if len(sort_rH_idx[shift_at[0]])<len(sort_rH_idx[shift_at[1]]): - changed from this line in 15032023
                transfer_type='R2L'
                tt_v=+1
                
        else:
            transfer_type='NA'
            tt_v=0
        
        
        ## update the index in the pH_idx (product H atom list) in line with the indexing in sort_rH_idx
        
        fix_pH_idx=pH_idx.copy()
    
        pls=[pH_idx[i] for i in aff_at]
        pls2=[]
        for ls in pls:
            pls2.append([i+tt_v for i in ls])
        
        for i,j in zip(aff_at,pls2):
            fix_pH_idx[i]=j
        
        #print(fix_pH_idx)
        
        ## find the atom involved in the tranfer 
        ## update the H atom index that is attached to the same non-H atom in the fix_pH_idx
        
        if transfer_type == 'L2R':
            tranfer_at=[i for i in sort_rH_idx[shift_at[0]] if i not in pH_idx[shift_at[0]]]
            if len(tranfer_at)>1:
                tranfer_at=[tranfer_at[0]]
    
            fix_pH_idx[shift_at[1]]=[i+tt_v for i in pH_idx[shift_at[1]] if len(pH_idx[shift_at[1]])>1  
                                    and i+tt_v not in naff_at_ls2]+tranfer_at
        
        elif transfer_type == 'R2L':
            tranfer_at=[i for i in sort_rH_idx[shift_at[1]] if i not in pH_idx[shift_at[1]]]
            if len(tranfer_at)>1:
                tranfer_at=[tranfer_at[0]]
            fix_pH_idx[shift_at[0]]=[i+tt_v for i in pH_idx[shift_at[0]] if len(pH_idx[shift_at[0]])>1  
                                    and i+tt_v not in naff_at_ls2]+tranfer_at


        
        ########
        ## the second filter 1. chk if the index of H atoms match in the two systems; correct if not 
       
        match_ls2=[]
        for i,j in zip(sort_rH_idx,fix_pH_idx):
            if i != j and len(i) != len(j):
                ## False,False -- H involved in the shift 
                match_ls2.append([False,False])
            elif i != j and len(i) == len(j):
                ## False,True -- not involved in the shift, but its indexing affected by the shift 
                match_ls2.append([False,True])
            else:
                ## True, True -- not involved in the shift and consistent index in sort_rH_idx and fix_pH_idx
                match_ls2.append([True,True])
        
        shift_at2=[i for i in range(0,len(match_ls)) if match_ls2[i] == [False,False]]
        aff_at2=[idx for idx in range(0,len(match_ls2)) if match_ls2[idx] == [False,True]]

        fix_pH_idx2=fix_pH_idx.copy()
        #print(fix_pH_idx2)
        
        ## insert the index of the shifting H into the correct sub list in fix_pH_idx2
        
        for idx in aff_at2:
            fix_pH_idx2[idx] = sort_rH_idx[idx]
        
        if transfer_type == 'L2R':

            fix_pH_idx2[shift_at2[0]] = [i for i in sort_rH_idx[shift_at2[0]] if i != tranfer_at[0]]
            fix_pH_idx2[shift_at2[1]] = sort_rH_idx[shift_at2[1]]+tranfer_at
         
        elif transfer_type == 'R2L':
            #print(tranfer_at)
            #print()
            fix_pH_idx2[shift_at2[0]] = sort_rH_idx[shift_at2[0]]+tranfer_at
            fix_pH_idx2[shift_at2[1]] =[i for i in sort_rH_idx[shift_at2[1]] if i !=tranfer_at[0]]
        
        
        
        ### finish sorting 
        ### constructing the dictionary now 
        
        pH_idx_ls=[i-1 for i in flatten(pH_idx)]
        fix_pH_idx_ls=flatten(fix_pH_idx2)
        
        
        rH_idx_ls= [i-1 for i in flatten(sort_rH_idx)]
        sort_rH_idx_ls=flatten(sort_rH_idx)
        
        
        self.pdict2=pdict.copy()
        self.rdict2=rdict.copy()
        
        for idx,midx in zip(pH_idx_ls,fix_pH_idx_ls):
            self.pdict2[idx]=midx
    
        for idx,midx in zip(rH_idx_ls,sort_rH_idx_ls):
            self.rdict2[idx]=midx
        
        #print(self.pdict2)
        #print(self.rdict2)
        
        
        self.fpr_dict={}
        self.frp_dict={}
        for i in range(0,len(self.rdict2)):    
            for j in range(0,len(self.rdict2)):
                #print(i,j)
                if self.rdict2.get(i)==self.pdict2.get(j):
                    #print(i,j)
                    self.fpr_dict[i]=j
                    self.frp_dict[j]=i
                    
 
        prod_H_at = [at for at in self.prod_H.GetAtoms()]

        for i in range(0,len(self.fpr_dict)):
            #prod_H_at[self.fpr_dict.get(i)].SetAtomMapNum(self.rdict2.get(i))
            prod_H_at[i].SetAtomMapNum(self.pdict2.get(i))

        
        react_H_at = [at for at in self.react_H.GetAtoms()]

        for i in range(0,len(self.fpr_dict)):
            react_H_at[i].SetAtomMapNum(self.rdict2.get(i))
        
        

    
    def get_info(self):
        
        
        ## note: to visualise the prediction heatmap, the input rxn smi string need to be
        ## the same as the input to the ML pipeline; otherwise, there might be inconsistence in the
        ## self.pdict2 and self.rdict2 dictionary 
        
        
        ## 
        
        react_symbol_ls=[at.GetSymbol() for at in self.react.GetAtoms()]
        prod_symbol_ls=[at.GetSymbol() for at in self.prod.GetAtoms()]
        
        self.smi_with_H=False
        
        if 'H' in react_symbol_ls and 'H' in prod_symbol_ls:
            self.smi_with_H=True
        
        
        if self.smi_with_H==False:
            self.process_rxn()
            
        else:
            self.react_H=self.react
            self.prod_H=self.prod
            
            react_at_mapidx=[at.GetAtomMapNum() for at in self.react_H.GetAtoms()]
            prod_at_mapidx=[at.GetAtomMapNum() for at in self.prod_H.GetAtoms()]
            
            self.rdict2={idx:react_at_mapidx[idx] for idx in range(0,len(react_at_mapidx))}
            self.pdict2={idx:prod_at_mapidx[idx] for idx in range(0,len(prod_at_mapidx))}
            
            self.fpr_dict={}
            self.frp_dict={}
            for i in range(0,len(self.rdict2)):    
                for j in range(0,len(self.rdict2)):
                    #print(i,j)
                    if self.rdict2.get(i)==self.pdict2.get(j):
                        #print(i,j)
                        self.fpr_dict[i]=j
                        self.frp_dict[j]=i
                        
        
################


def mapping_competitive_pathways(reaction_df, v):
    
    sub_rxn_df=reaction_df[reaction_df['code']==v]
    rxn_ls=sub_rxn_df['reaction'].tolist()
    idx_ls=sub_rxn_df['idx'].tolist()
    code_ls=sub_rxn_df['code'].tolist()
    #ea_ls=sub_rxn_df['ea'].tolist()
    
    code_value=code_ls[0]
    
    rxn=rxn_ls[0]
    
    reaction=AtIdx([rxn])
    reaction.get_info()
    
    react_smi=Chem.MolToSmiles(reaction.react_H)
    prod_smi=Chem.MolToSmiles(reaction.prod_H)
    
    new_rxn_ls=[react_smi+'>>'+prod_smi]
    
    error_ls=[False]
    
    
    for idx in range(1,len(rxn_ls)):
        error_tag=False
        rxn2=rxn_ls[idx]
        
        reaction2=AtIdx([rxn2])
        reaction2.get_info()
        
        mapping=list(reaction.react_H.GetSubstructMatch(reaction2.react_H))

        reidx_dict={}
        
        for atom,i in zip(reaction2.react_H.GetAtoms(),mapping):
            original=atom.GetAtomMapNum()
            map_idx=reaction.react_H.GetAtomWithIdx(i).GetAtomMapNum()
            atom.SetAtomMapNum(map_idx)
            reidx_dict[original]=map_idx
            
        for atom in reaction2.prod_H.GetAtoms():
            ori=atom.GetAtomMapNum()
            atom.SetAtomMapNum(reidx_dict.get(ori))
            
        react_smi2=Chem.MolToSmiles(reaction2.react_H)
        prod_smi2=Chem.MolToSmiles(reaction2.prod_H)
        
        if react_smi2 != react_smi:
            error_tag=True
            print('mismatch reactant smi string, code:'+ str(code_value)+' idx:' 
                  + str(idx))
            
        new_rxn_ls.append(react_smi2+'>>'+prod_smi2)
        error_ls.append(error_tag)
    
    #new_rxn_df=pd.DataFrame({'idx':idx_ls,'code':code_ls, 'reaction':new_rxn_ls,'ea':ea_ls,
    #                        'error':error_ls})

    new_rxn_df=pd.DataFrame({'idx':idx_ls,'code':code_ls, 'reaction':new_rxn_ls,'error':error_ls})
    
    return new_rxn_df
    
    
    
def mapping_competitive_pathways_all(csv):
    
    reaction_df = pd.read_csv(csv)
    
    unique_code = reaction_df['code'].unique().tolist()
    sub_rxn_df_ls=[]
    for v in unique_code:
        #try:
        result=mapping_competitive_pathways(reaction_df, v)
        sub_rxn_df_ls.append(result)
        #except:
        #    print(v)


    rxn_df = pd.concat(sub_rxn_df_ls)

    rxn_df['code'] = parah_code(rxn_df['code'].tolist())
    rxn_df['idx'] = [i for i in range(0,len(rxn_df))]
    
    return rxn_df 


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
            
            
            
            
            