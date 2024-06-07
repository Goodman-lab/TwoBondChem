#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon July 03 14:05:13 2023

@author: chingchinglam
"""

from datetime import date
import sys
import os
sys.path.insert(1, os.getcwd()+'/scripts/')


import pandas as pd
import model_training_v3 as mt
import random
            
def give_random_sample(filename):
    data = pd.read_csv(filename)

    random_numbers = []
    while len(random_numbers) < 100:
        va=random.randint(0, data['code'].max())
        if va not in random_numbers:
            random_numbers.append(va)


    random_code_ls=sorted(random_numbers)
    #print(len(random_code_ls))
    sample_df=data.loc[data['code'].isin(random_code_ls)].copy()

    sample_df['code'] = parah_code(sample_df['code'].tolist())
    sample_df['idx'] = [i for i in range(0,len(sample_df))]

    return sample_df

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


today = date.today()
date_str=today.strftime("%d%m%Y")
this_python_script_name= sys.argv[0].split('/')[-1][:-3]


###################

filename_ls = ['diels_alder_data_v7_19052024.csv','cyclo_data_v2_13012024.csv']
ram_data_ls = [give_random_sample(f) for f in filename_ls]
ele_steps = pd.concat(ram_data_ls)

ele_steps['code'] = parah_code(ele_steps['code'].tolist())
ele_steps['idx'] = [i for i in range(0,len(ele_steps))]


test=mt.run_ML(ele_steps)
test.generate_data()
print('done: test.generate_data()')

test.training(save_model = True, file_name='RF_model_Bmix_')


###################

filename_ls = ['fav_RGD1_13012024_v2.csv','green_04012024.csv']
ram_data_ls = [give_random_sample(f) for f in filename_ls]
ele_steps = pd.concat(ram_data_ls)

ele_steps['code'] = parah_code(ele_steps['code'].tolist())
ele_steps['idx'] = [i for i in range(0,len(ele_steps))]


test=mt.run_ML(ele_steps)
test.generate_data()
print('done: test.generate_data()')

test.training(save_model = True, file_name='RF_model_Cmix_')


###################

yr1_data=pd.read_csv('first-year_data_all_27122023.csv')

filename_ls = ['diels_alder_data_v7_19052024.csv','cyclo_data_v2_13012024.csv',
'fav_RGD1_13012024_v2.csv','green_04012024.csv']



ram_data_ls = [give_random_sample(f) for f in filename_ls]
ele_steps = pd.concat([yr1_data]+ram_data_ls)

ele_steps['code'] = parah_code(ele_steps['code'].tolist())
ele_steps['idx'] = [i for i in range(0,len(ele_steps))]


test=mt.run_ML(ele_steps)
test.generate_data()
print('done: test.generate_data()')

test.training(save_model = True, file_name='RF_model_all_')




