#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 09:55:45 2024

@author: chingchinglam
"""


import pandas as pd 
import json 
from sklearn.metrics import accuracy_score,precision_score,recall_score
import os 


def get_label_ls(data_df):
    actual_label_ls=[]
    #print(data_df['actual_label'])

    for ls in data_df['actual_label'].tolist():
        actual_label_ls.extend(json.loads(ls))
        
    pred_label_ls=[]

    for ls in data_df['pred_label'].tolist():
        pred_label_ls.extend(json.loads(ls))
        
    return actual_label_ls,pred_label_ls

def get_label_ls_v2(data_df):
    actual_label_ls=[]

    for ls in data_df['actual_label'].tolist():
        actual_label_ls.append(json.loads(ls))
        
    pred_label_ls=[]

    for ls in data_df['pred_label'].tolist():
        pred_label_ls.append(json.loads(ls))
        
    return actual_label_ls,pred_label_ls



def get_metrice(pred_ls,test_ls):
    
    accuracy = accuracy_score(test_ls,pred_ls)
    recall = recall_score(test_ls,pred_ls)
    precision = precision_score(test_ls,pred_ls)

    label_ls=[]

    for p,t in zip(pred_ls,test_ls):
        if p == 1 and t == 1:
            label_ls.append('tp')
        elif p == 0 and t == 0:
            label_ls.append('tn')
        elif p == 1 and t == 0:
            label_ls.append('fp')
        elif p == 0 and t == 1:
            label_ls.append('fn')
        
    r_tp = len([i for i in label_ls if i == 'tp']) / len(label_ls)
    r_tn = len([i for i in label_ls if i == 'tn']) / len(label_ls)
    r_fp = len([i for i in label_ls if i == 'fp']) / len(label_ls)
    r_fn = len([i for i in label_ls if i == 'fn']) / len(label_ls)
    
    no_fp = len([i for i in label_ls if i == 'fp'])
    no_fn = len([i for i in label_ls if i == 'fn'])
    
    
    True_label_pred=len([i for i in pred_ls if i >0])/len(pred_ls)
    True_label_act=len([i for i in test_ls if i >0])/len(test_ls)

    result_dict={'accuracy':accuracy,'recall':recall, 'precision':precision,'tp':r_tp,'tn':r_tn, 
                'fp':r_fp, 'fn':r_fn,'n_fp':no_fp, 'n_fn':no_fn,'postive_label_pred':True_label_pred,
                'postive_label_act':True_label_act}
    
    return result_dict





def comp_evaluation(opt_files_raw):
    
    opt_files_name = ['_'.join(f.split('_')[0:-1]) for f in opt_files_raw]

    accur_ls=[]
    recall_ls=[]
    precision_ls=[]
    
    fp_ls=[]
    fn_ls=[]
    
    nfn_ls=[]
    no_fault_ls=[]
    no_fault1_ls=[]
    
    positive_label_ls=[]
    
    for csv in opt_files_raw:
        data = pd.read_csv(csv)
        act_lab,pred_lab=get_label_ls(data)
        
        metrice_dict = get_metrice(pred_lab,act_lab)
    
        accur_ls.append(metrice_dict.get('accuracy'))
        precision_ls.append(metrice_dict.get('precision'))
        fp_ls.append(metrice_dict.get('fp'))
        fn_ls.append(metrice_dict.get('fn'))
        positive_label_ls.append(metrice_dict.get('postive_label_act'))
        recall_ls.append(metrice_dict.get('recall'))
    
    
        rx_act_lab,rx_pred_lab=get_label_ls_v2(data)
    
        xfn_ls=[]
        xfp_ls=[]
        xallfault_ls=[]
        
        for act_ls,pred_ls in zip(rx_act_lab,rx_pred_lab):
    
            metrice_dict = get_metrice(pred_ls,act_ls)
    
            fn=metrice_dict.get('n_fn')
            fp=metrice_dict.get('n_fp')
            xfn_ls.append(fn)
            xfp_ls.append(fp)
            xallfault_ls.append(fn+fp)
            
    
        nfn_ls.append(len([i for i in xfn_ls if i == 0])/len(xfn_ls))
    
        no_fault_ls.append(len([i for i in xallfault_ls if i == 0])/len(xallfault_ls))
        no_fault1_ls.append(len([i for i in xallfault_ls if i <= 1])/len(xallfault_ls))
    
        
        
        
    result=pd.DataFrame({'file':opt_files_name, 'accuracy':accur_ls,  'recall':recall_ls,'precision':precision_ls,
                         'no_fault':no_fault_ls,'no_fault1':no_fault1_ls, 'n_fn':nfn_ls})
    
    
    return result 



#opt_files_raw = sorted([f for f in os.listdir('./') if f.endswith('.csv')])





