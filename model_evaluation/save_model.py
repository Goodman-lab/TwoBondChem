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
            
                

today = date.today()
date_str=today.strftime("%d%m%Y")
this_python_script_name= sys.argv[0].split('/')[-1][:-3]


ele_steps=pd.read_csv('cyclo_data_v2_13012024.csv').head(100)

test=mt.run_ML(ele_steps)
test.generate_data()
print('done: test.generate_data()')

test.training(save_model = True, file_name='RF_model_cyclo_')

ele_steps=pd.read_csv('diels_alder_data_v7_19052024.csv').head(100)

test=mt.run_ML(ele_steps)
test.generate_data()
print('done: test.generate_data()')

test.training(save_model = True, file_name='RF_model_da_')

ele_steps=pd.read_csv('fav_RGD1_13012024_v2.csv').head(300)

test=mt.run_ML(ele_steps)
test.generate_data()
print('done: test.generate_data()')

test.training(save_model = True, file_name='RF_model_rgd_')



ele_steps=pd.read_csv('green_04012024.csv').head(100)

test=mt.run_ML(ele_steps)
test.generate_data()
print('done: test.generate_data()')

test.training(save_model = True, file_name='RF_model_green_')


ele_steps=pd.read_csv('first-year_data_all_27122023.csv')

test=mt.run_ML(ele_steps)
test.generate_data()
print('done: test.generate_data()')

test.training(save_model = True, file_name='RF_model_yr1_')


