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

import pandas as pd
sys.path.insert(1, os.getcwd()+'/scripts/')


import model_training_v3 as mt


today = date.today()
date_str=today.strftime("%d%m%Y")
this_python_script_name= sys.argv[0].split('/')[-1][:-3]

#test_df = test_df1[test_df1['code']<100]


def Onepot(training_df, test_df, model ='RF', descriptor_arr = '2-bond+'):
    
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
    test.generate_data(method='m3',descriptor_arr =descriptor_arr)
    print('done: test.generate_data()')
    #test.status_df.to_csv('status.csv')
    
    training_code_ls=[test.rev_pass_idx_dict.get(i) for i in ele_steps[ele_steps['type']=='train']['code'].unique().tolist() 
                      if i in test.pass_idx_ls ]
    test_code_ls=[test.rev_pass_idx_dict.get(i) for i in ele_steps[ele_steps['type']=='test']['code'].unique().tolist() 
                  if i in test.pass_idx_ls ]
    
    how_ls=[training_code_ls,test_code_ls]
    
    
    test.split_data(by='rxn',incl='r',test_size=0.2, how=how_ls)
    
    #print('idx: ', str(looping_ls[idx][1]))
    
    test.training(model=model)
    test.evaluation()
    
    
    original_code_ls = [i-len_training_code-1 for i in test.label_df['rxn_idx'].tolist()]
    #original_idx_ls = [i-len_training_idx-1 for i in test.label_df['idx'].tolist()]
    
    label_df2=test.label_df[['pred_label','actual_label']].copy()
    
    label_df2['rxn_idx']=original_code_ls
    label_df2['idx']=original_code_ls
    
    return label_df2
    
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

def give_random_sample_yr1(filename):
    data = pd.read_csv(filename)

    random_numbers = []
    while len(random_numbers) < 78:
        va=random.randint(0, data['code'].max())
        if va not in random_numbers:
            random_numbers.append(va)

    

    random_code_ls=sorted(random_numbers)

    full_code_ls = data['code'].unique().tolist()
    test_code_ls = [i for i in full_code_ls if i not in random_code_ls]


    train_sample_df=data.loc[data['code'].isin(random_code_ls)].copy()

    train_sample_df['code'] = parah_code(train_sample_df['code'].tolist())
    train_sample_df['idx'] = [i for i in range(0,len(train_sample_df))]

    test_sample_df=data.loc[data['code'].isin(test_code_ls)].copy()

    test_sample_df['code'] = parah_code(test_sample_df['code'].tolist())
    test_sample_df['idx'] = [i for i in range(0,len(test_sample_df))]

    return train_sample_df, test_sample_df


def give_random_sample(filename):
    data = pd.read_csv(filename)

    random_numbers = []
    while len(random_numbers) < 200:
        va=random.randint(0, data['code'].max())
        if va not in random_numbers:
            random_numbers.append(va)


    random_code_ls=sorted(random_numbers)
    #print(len(random_code_ls))
    sample_df=data.loc[data['code'].isin(random_code_ls)].copy()

    sample_df['code'] = parah_code(sample_df['code'].tolist())
    sample_df['idx'] = [i for i in range(0,len(sample_df))]

    return sample_df


filename = './dataset/first-year_data_all_27122023.csv'


filename_ls = ['./dataset/cyclo_data_v2_16072024.csv', './dataset/diels_alder_data_v7_19052024.csv', 
                './dataset/fav_RGD1_13012024_v2_nodu.csv', './dataset/green_04012024.csv']

keyword = 'all_2b+'

model_ls = []

test_df_ls=[]
train_df_ls=[]


#'cyclo','da',,'RGD1','green'

for i in range(0,10):

    model_ls.append('_'+str(i))

    training_df_yr1, test_df_yr1 =give_random_sample_yr1(filename)

    training_df_raw_ls=[training_df_yr1]
    test_df_raw_ls=[test_df_yr1]

    for fil in filename_ls:
        sample_df1=give_random_sample(fil)
        test_df1 = sample_df1[sample_df1['code']<100]
        training_df1 = sample_df1[sample_df1['code']>=100]
        training_df_raw_ls.append(training_df1)
        test_df_raw_ls.append(test_df1)

    comb_training_df=pd.concat(training_df_raw_ls)
    comb_testing_df=pd.concat(test_df_raw_ls)

    test_df=comb_testing_df.copy()
    train_df=comb_training_df.copy()

    test_df['code'] = parah_code(comb_testing_df['code'].tolist())
    test_df['idx'] = [i for i in range(0,len(comb_testing_df))]

    train_df['code'] = parah_code(comb_training_df['code'].tolist())
    train_df['idx'] = [i for i in range(0,len(comb_training_df))]

    test_df_ls.append(test_df)
    train_df_ls.append(train_df)



label_df_ls=[]
for test,train,mod in zip(test_df_ls,train_df_ls,model_ls):
    result_df=Onepot(train,test, model = 'RF',descriptor_arr = '2-bond+')
    result_df.to_csv(keyword +mod+'_'+date_str+'.csv')
    label_df_ls.append(result_df)


keyword = 'all_2b'

label_df_ls=[]
for test,train,mod in zip(test_df_ls,train_df_ls,model_ls):
    result_df=Onepot(train,test, model = 'RF',descriptor_arr = '2-bond')
    result_df.to_csv(keyword +mod+'_'+date_str+'.csv')
    label_df_ls.append(result_df)



keyword = 'all_1b'

label_df_ls=[]
for test,train,mod in zip(test_df_ls,train_df_ls,model_ls):
    result_df=Onepot(train,test, model = 'RF',descriptor_arr = '1-bond')
    result_df.to_csv(keyword +mod+'_'+date_str+'.csv')
    label_df_ls.append(result_df)







