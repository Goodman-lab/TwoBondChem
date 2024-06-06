#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 11:29:11 2022

@author: chingchinglam
"""


import pandas as pd
import numpy as np
import scipy.stats as stats
import math
from collections import defaultdict
import os 
from datetime import date
import sys
import json

def list_duplicates(seq):

    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)  
    
    return [[key,len(locs)] for key,locs in tally.items() if len(locs)>0]

def percentile(data, perc: int):
    size = len(data)
    return sorted(data)[int(math.ceil((size * perc) / 100)) - 1]


class class_bond:
    def __init__(self, result):
        self.unique_bond=result['sym'].unique().tolist()
        self.sub_df_ls=[result[result['sym']==b] for b in self.unique_bond]
        
    def get_boundary(self, bondset):
        bondlength_ls=np.array(bondset['length'])

        b_max=max(bondlength_ls)+0.1
        b_min=min(bondlength_ls)-0.1
        #b_max = percentile(bondlength_ls,98)
        #b_min = percentile(bondlength_ls,1)

        sel_num=round((b_max-b_min)*2000)

        frame_ls= np.linspace(b_min, b_max, num=sel_num)

        kde = stats.gaussian_kde(bondlength_ls)
        kde.set_bandwidth(bw_method=kde.factor / 3)
        #kde.set_bandwidth(bw_method='silverman')

        y1 = kde(frame_ls)

        gy1=np.gradient(y1)

        s = np.sign(gy1)


        boundary_ls=[0]
        s_ls=s.tolist()
        init=s_ls[0]

        for idx in range(1,len(s_ls)):
            
            if s_ls[idx] != init:
                boundary_ls.append(idx)
                init = s_ls[idx]
                
        
        boundary_ls1 =   [frame_ls[i] for i in boundary_ls] + [frame_ls[-1]]
        
        class_data=pd.cut(bondset.length, boundary_ls1, labels=[str(i) for i in list(range(1,len(boundary_ls1)))],
                 duplicates='drop').tolist()
        
        
        sort_ls=list_duplicates(class_data)
        sort_ls2 = [i for i in sort_ls if i[1] > 10]


        select_idx = sorted([int(i[0])-1 for i in sort_ls2])

        if len(sort_ls)-1 != select_idx[-1]:  
            select_idx2 = select_idx+[select_idx[-1]+1]
        else:
            select_idx2 = select_idx

        boundary_ls2 = [boundary_ls1[0]]+[boundary_ls1[idx] for idx in select_idx2]+[boundary_ls1[-1]]

        boundary_select = [boundary_ls[i] for i in select_idx2]
        
        
        max_or_min='max'        
        max_or_min_ls=[]
        for i in boundary_select:
            if gy1[i-2] > 0 and gy1[i+2] <0:
                max_or_min_ls.append('max')
            else:
                max_or_min_ls.append('min')
    
        min_idx_ls=[idx for idx in range(0,len(max_or_min_ls)) if max_or_min_ls[idx] == 'min']
    
        select_idx3=[select_idx2[idx] for idx in min_idx_ls]

        self.boundary_ls3 = [boundary_ls1[0]]+[boundary_ls1[idx] for idx in select_idx3]+[boundary_ls1[-1]]
        
        class_data2=pd.cut(bondset.length, self.boundary_ls3, labels=[str(i) for i in list(range(1,len(set(self.boundary_ls3))))],
                           duplicates='drop').tolist()
        self.sort_ls_fin=list_duplicates(class_data2)
    
        
    def perform(self):
        
        boundary_info_ls=[]
        sortls_ls=[]
        bond_ls=[]
        for df,b in zip(self.sub_df_ls,self.unique_bond):
            try:
                self.get_boundary(df)
                boundary_info_ls.append(self.boundary_ls3 )
                sortls_ls.append(self.sort_ls_fin)
                bond_ls.append(b)
            except:
                print(b)
            
            
            
        self.result_df=pd.DataFrame({'bond':bond_ls,
                                'boundary':boundary_info_ls,
                                'sort':sortls_ls})



def get_bond_classification(result):

    boundary_ls=[list(set(json.loads(i))) for i in result['boundary'].tolist()]

    sort_ls =  result['sort'].tolist()
    sort_ls2 = [i.split("'") for i in sort_ls] 
    sort_ls3 = [''.join(i) for i in sort_ls2 ]
    sort_ls4 = [sorted(json.loads(i)) for i in sort_ls3]

    bond_ls=result['bond']

    sort_ls5=[]
    bond_ls2=[]
    boundary_ls2=[]
    for sls,b,bd in zip(sort_ls4,bond_ls,boundary_ls):
        if len(sls) >1:
            sort_ls5.append(sls)
            bond_ls2.append(b)
            boundary_ls2.append(bd)

            
    result2 = pd.DataFrame({'bond':bond_ls2,'boundary':boundary_ls2,'sort':sort_ls5})

    xbond_ls=['NN','NS','BrC']
    fresult=pd.concat([result2[result2['bond']==i] for i in result2['bond'].unique() if i not in xbond_ls])

    return fresult





path = os.getcwd()+'/data/'
today = date.today()
date_str=today.strftime("%d%m%Y")
this_python_script_name= sys.argv[0].split('/')[-1][:-3]

result1=pd.read_csv('bondlength_fromsmi_v2_.csv')

test2=class_bond(result1)
test2.perform()
result2=get_bond_classification(test2.result_df)

result2.to_csv('bond_classification_'+date_str+'.csv')

print('Done: '+this_python_script_name+'_7_'+date_str+'.csv')







