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

def give_random_sample(filename,sample_size=200):
    data = pd.read_csv(filename)

    random_numbers = []
    while len(random_numbers) < sample_size:
        va=random.randint(0, data['code'].max())
        if va not in random_numbers:
            random_numbers.append(va)


    random_code_ls=sorted(random_numbers)
    #print(len(random_code_ls))
    sample_df=data.loc[data['code'].isin(random_code_ls)].copy()

    sample_df['code'] = parah_code(sample_df['code'].tolist())
    sample_df['idx'] = [i for i in range(0,len(sample_df))]

    return sample_df



filename = './dataset/fav_RGD1_13012024_v2.csv'
keyword = 'rgd_size'

## uncomment the below to conduct tests on these dataset 

#filename = './dataset/diels_alder_data_v7_19052024.csv'
#keyword = 'da_size'

#filename = './dataset/cyclo_data_v2_16072024.csv'
#keyword = 'cycl_size'


model_ls = []

test_df_ls=[]
train_df_ls=[]


test_size_ls=[100,200,500,1000,1500,2000]


#'cyclo','da',,'RGD1','green'

for j in range(0,5):

    for i in test_size_ls:

        model_ls.append('_'+str(i)+'_'+str(j))

        sample_df=give_random_sample(filename,sample_size=i+100)

        test_df = sample_df[sample_df['code']<i]
        training_df = sample_df[(sample_df['code']>=i)]

        ### uncomment the below if you wish to vary the training data set instead
        #training_df = sample_df[sample_df['code']<i]
        #test_df = sample_df[sample_df['code']<i]

        test_df_ls.append(test_df)
        train_df_ls.append(training_df)



label_df_ls=[]
for test,train,mod in zip(test_df_ls,train_df_ls,model_ls):
    result_df=Onepot(train,test, model = 'RF',descriptor_arr = '2-bond+')
    result_df.to_csv(keyword +mod+'_'+date_str+'.csv')
    label_df_ls.append(result_df)








