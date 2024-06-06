#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 11:19:31 2023

@author: chingchinglam
"""

from datetime import date
import sys
import os

import pandas as pd
sys.path.insert(1, os.getcwd()+'/scripts/')


import model_training_v3 as mt


today = date.today()
date_str=today.strftime("%d%m%Y")
this_python_script_name= sys.argv[0].split('/')[-1][:-3]

#####################

ele_stepsx = pd.read_csv('cyclo_data_v2_04012024.csv')
ele_steps=ele_stepsx[ele_stepsx['code']<=100]



test=mt.run_ML(ele_steps)
test.generate_data(method='m3')
print('done: test.generate_data()')


looping_ls = mt.get_LoopingLs(test.pass_sub_df,v=1)


result_df_ls=[]
label_df_ls=[]
#shap_value_df_ls=[]

for idx in range(0, len(looping_ls)):
    test.split_data(by='rxn',incl='r',test_size=0.2, how=looping_ls[idx])
    
    #print('idx: ', str(looping_ls[idx][1]))
    
    test.training(model='classifier')
    test.evaluation()
    #test.XAL_test()
    
    test.descriptor_df['rxn_idx']=idx
    test.descriptor_df['rxn_idx']=idx
    
    result_df_ls.append(test.descriptor_df)
    label_df_ls.append(test.label_df)
    #shap_value_df_ls.append(test.shap_df)
    
    
    

result_df=pd.concat(result_df_ls)
label_df=pd.concat(label_df_ls)

label_df.to_csv(this_python_script_name+'_'+date_str+'_label_cylo.csv')


#####################


ele_stepsx = pd.read_csv('fav_RGD1_04012024.csv')
ele_steps=ele_stepsx[ele_stepsx['code']<=100]


test=mt.run_ML(ele_steps)
test.generate_data(method='m3')
print('done: test.generate_data()')

looping_ls = mt.get_LoopingLs(test.pass_sub_df,v=1)


result_df_ls=[]
label_df_ls=[]
#shap_value_df_ls=[]

for idx in range(0, len(looping_ls)):
    test.split_data(by='rxn',incl='r',test_size=0.2, how=looping_ls[idx])
    
    #print('idx: ', str(looping_ls[idx][1]))
    
    test.training(model='classifier')
    test.evaluation()
    #test.XAL_test()
    
    test.descriptor_df['rxn_idx']=idx
    test.descriptor_df['rxn_idx']=idx
    
    result_df_ls.append(test.descriptor_df)
    label_df_ls.append(test.label_df)
    #shap_value_df_ls.append(test.shap_df)
    
    

result_df=pd.concat(result_df_ls)
label_df=pd.concat(label_df_ls)


label_df.to_csv(this_python_script_name+'_'+date_str+'_label_rgd.csv')

#####################

ele_stepsx = pd.read_csv('diels_alder_data_v7_19052024.csv')
ele_steps=ele_stepsx[ele_stepsx['code']<=100]
#rxns_ls2 = ele_steps['reaction'].tolist()

test=mt.run_ML(ele_steps)
test.generate_data(method='m3')
print('done: test.generate_data()')

looping_ls = mt.get_LoopingLs(test.pass_sub_df,v=1)


result_df_ls=[]
label_df_ls=[]
#shap_value_df_ls=[]

for idx in range(0, len(looping_ls)):
    test.split_data(by='rxn',incl='r',test_size=0.2, how=looping_ls[idx])
    
    #print('idx: ', str(looping_ls[idx][1]))
    
    test.training(model='classifier')
    test.evaluation()
    #test.XAL_test()
    
    test.descriptor_df['rxn_idx']=idx
    test.descriptor_df['rxn_idx']=idx
    
    result_df_ls.append(test.descriptor_df)
    label_df_ls.append(test.label_df)
    #shap_value_df_ls.append(test.shap_df)
    
    

result_df=pd.concat(result_df_ls)
label_df=pd.concat(label_df_ls)
#shap_value_df=pd.concat(shap_value_df_ls)


label_df.to_csv(this_python_script_name+'_'+date_str+'_label_da.csv')


#####################

#####################

ele_stepsx = pd.read_csv('da_08012024_vCorr.csv')
ele_steps=ele_stepsx[ele_stepsx['code']<=100]
#rxns_ls2 = ele_steps['reaction'].tolist()

test=mt.run_ML(ele_steps)
test.generate_data(method='m3')
print('done: test.generate_data()')

looping_ls = mt.get_LoopingLs(test.pass_sub_df,v=1)


result_df_ls=[]
label_df_ls=[]
#shap_value_df_ls=[]

for idx in range(0, len(looping_ls)):
    test.split_data(by='rxn',incl='r',test_size=0.2, how=looping_ls[idx])
    
    #print('idx: ', str(looping_ls[idx][1]))
    
    test.training(model='classifier')
    test.evaluation()
    #test.XAL_test()
    
    test.descriptor_df['rxn_idx']=idx
    test.descriptor_df['rxn_idx']=idx
    
    result_df_ls.append(test.descriptor_df)
    label_df_ls.append(test.label_df)
    #shap_value_df_ls.append(test.shap_df)
    
    

result_df=pd.concat(result_df_ls)
label_df=pd.concat(label_df_ls)
#shap_value_df=pd.concat(shap_value_df_ls)


label_df.to_csv(this_python_script_name+'_'+date_str+'_label_da_vCorr.csv')


#####################

ele_stepsx = pd.read_csv('green_04012024.csv')
ele_steps=ele_stepsx[ele_stepsx['code']<=100]
#rxns_ls2 = ele_steps['reaction'].tolist()



test=mt.run_ML(ele_steps)
test.generate_data(method='m3')
print('done: test.generate_data()')


looping_ls = mt.get_LoopingLs(test.pass_sub_df,v=1)


result_df_ls=[]
label_df_ls=[]
#shap_value_df_ls=[]

for idx in range(0, len(looping_ls)):
    test.split_data(by='rxn',incl='r',test_size=0.2, how=looping_ls[idx])
    
    #print('idx: ', str(looping_ls[idx][1]))
    
    test.training(model='classifier')
    test.evaluation()
    #test.XAL_test()
    
    test.descriptor_df['rxn_idx']=idx
    test.descriptor_df['rxn_idx']=idx
    
    result_df_ls.append(test.descriptor_df)
    label_df_ls.append(test.label_df)
    #shap_value_df_ls.append(test.shap_df)
    
    

result_df=pd.concat(result_df_ls)
label_df=pd.concat(label_df_ls)

label_df.to_csv(this_python_script_name+'_'+date_str+'_label_green.csv')

#####################


ele_steps=pd.read_csv('first-year_data_all_27122023.csv')


test=mt.run_ML(ele_steps)
test.generate_data(method='m3')
print('done: test.generate_data()')


looping_ls = mt.get_LoopingLs(test.pass_sub_df,v=1)


result_df_ls=[]
label_df_ls=[]
#shap_value_df_ls=[]

for idx in range(0, len(looping_ls)):
    test.split_data(by='rxn',incl='r',test_size=0.2, how=looping_ls[idx])
    
    #print('idx: ', str(looping_ls[idx][1]))
    
    test.training(model='classifier')
    test.evaluation()
    #test.XAL_test()
    
    test.descriptor_df['rxn_idx']=idx
    test.descriptor_df['rxn_idx']=idx
    
    result_df_ls.append(test.descriptor_df)
    label_df_ls.append(test.label_df)
    #shap_value_df_ls.append(test.shap_df)
    
    

result_df=pd.concat(result_df_ls)
label_df=pd.concat(label_df_ls)

label_df.to_csv(this_python_script_name+'_'+date_str+'_label_yr1.csv')

