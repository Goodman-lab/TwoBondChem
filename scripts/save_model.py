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
import model_training_v2 as mt
            
                

today = date.today()
date_str=today.strftime("%d%m%Y")
this_python_script_name= sys.argv[0].split('/')[-1][:-3]


ele_steps=pd.read_csv('first-year_data_all_27122023.csv')

test=mt.run_ML(ele_steps)
test.generate_data()
print('done: test.generate_data()')

test.training(save_model = True, file_name='RF_model_yr1_')






