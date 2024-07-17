#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 19:53:18 2023

@author: chingchinglam
"""

from datetime import date
import random
import sys
import os
import ast

import pandas as pd
sys.path.insert(1, os.getcwd()+'/scripts/')


import model_training_v3 as mt


today = date.today()
date_str=today.strftime("%d%m%Y")
this_python_script_name= sys.argv[0].split('/')[-1][:-3]

#test_df = test_df1[test_df1['code']<100]


def Onepot(training_df,test_df, model ='RF',
    hyperpara={'n_estimators':100,'max_features':'sqrt','criterion':'gini'}):
    
    training_df['type'] = 'train'
    test_df['type'] = 'test'
    
    
    len_training_code = max(training_df['code'].tolist())
    len_training_idx = max(training_df['idx'].tolist())
    
    test_code_ls = [i+len_training_code+1 for i in test_df['code'].tolist()]
    test_idx_ls = [i+len_training_idx+1 for i in test_df['idx'].tolist()]
    
    test_df2=test_df[['type','reaction']].copy()
    
    test_df2['code'] = test_code_ls
    test_df2['idx'] = test_idx_ls
    
    ele_steps=pd.concat([training_df,test_df2])

    ele_steps.to_csv('ele_steps.csv')
    
    test=mt.run_ML(ele_steps)
    test.generate_data(method='m3')
    print('done: test.generate_data()')
    #test.status_df.to_csv('status.csv')
    
    training_code_ls=[test.rev_pass_idx_dict.get(i) for i in ele_steps[ele_steps['type']=='train']['code'].unique().tolist() 
                      if i in test.pass_idx_ls ]
    test_code_ls=[test.rev_pass_idx_dict.get(i) for i in ele_steps[ele_steps['type']=='test']['code'].unique().tolist() 
                  if i in test.pass_idx_ls ]
    
    how_ls=[training_code_ls,test_code_ls]
    
    
    test.split_data(by='rxn',incl='r',test_size=0.2, how=how_ls)
    
    #print('idx: ', str(looping_ls[idx][1]))
    
    test.training(model=model,hyperpara=hyperpara)
    test.evaluation()
    
    
    original_code_ls = [i-len_training_code-1 for i in test.label_df['rxn_idx'].tolist()]
    #original_idx_ls = [i-len_training_idx-1 for i in test.label_df['idx'].tolist()]
    
    label_df2=test.label_df[['pred_label','actual_label']].copy()
    
    label_df2['rxn_idx']=original_code_ls
    label_df2['idx']=original_code_ls
    
    return label_df2


para_df = pd.read_csv('para_df.csv')

dict_ls=[]
para_str_ls=[]
for dit in para_df['params'].tolist():
    
    sdict=ast.literal_eval(dit)
    dict_ls.append(sdict)
    
    para_str_ls.append('_'.join([str(i) for i in list(sdict.values())]))

filename = './dataset/first-year_data_all_27122023.csv'
keyword = 'yr1'

training_df1 = pd.read_csv(filename)

test_df = training_df1[training_df1['code']>=78]
training_df = training_df1[training_df1['code']<78]




label_df_ls=[]
for para_dict, para_str in zip(dict_ls,para_str_ls):
    result_df=Onepot(training_df,test_df, model = 'RF',hyperpara=para_dict)
    result_df.to_csv(keyword +'_'+para_str+'_'+date_str+'.csv')
    label_df_ls.append(result_df)


#################

para_df = pd.read_csv('para_df.csv')

dict_ls=[]
para_str_ls=[]
for dit in para_df['params'].tolist():
    
    sdict=ast.literal_eval(dit)
    dict_ls.append(sdict)
    
    para_str_ls.append('_'.join([str(i) for i in list(sdict.values())]))

filename = './dataset/fav_RGD1_13012024_v2.csv'
keyword = 'rgd'

training_df1 = pd.read_csv(filename)

test_df = training_df1[training_df1['code']<100]
training_df = training_df1[(training_df1['code']>=100)&(training_df1['code']<200)]




label_df_ls=[]
for para_dict, para_str in zip(dict_ls,para_str_ls):
    result_df=Onepot(training_df,test_df, model = 'RF',hyperpara=para_dict)
    result_df.to_csv(keyword +'_'+para_str+'_'+date_str+'.csv')
    label_df_ls.append(result_df)


#################


para_df = pd.read_csv('para_df.csv')

dict_ls=[]
para_str_ls=[]
for dit in para_df['params'].tolist():
    
    sdict=ast.literal_eval(dit)
    dict_ls.append(sdict)
    
    para_str_ls.append('_'.join([str(i) for i in list(sdict.values())]))

filename = './dataset/cyclo_data_v2_16072024.csv'
keyword = 'cy'

training_df1 = pd.read_csv(filename)

test_df = training_df1[training_df1['code']<100]
training_df = training_df1[(training_df1['code']>=100)&(training_df1['code']<200)]




label_df_ls=[]
for para_dict, para_str in zip(dict_ls,para_str_ls):
    result_df=Onepot(training_df,test_df, model = 'RF',hyperpara=para_dict)
    result_df.to_csv(keyword +'_'+para_str+'_'+date_str+'.csv')
    label_df_ls.append(result_df)

