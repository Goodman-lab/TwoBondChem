#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon July 03 14:05:13 2023

@author: chingchinglam

## derive from test_looping_v5 
## renamed to call model_training.py

"""



from datetime import date
import sys
import os
import pickle
sys.path.insert(1, os.getcwd()+'/scripts/')

#from rdkit import Chem
import pandas as pd
import training_prepare_v5 as tp
import numpy as np
#from scipy.stats import spearmanr
from collections import defaultdict
import json 
from sklearn.metrics import accuracy_score,precision_score,recall_score

from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
#from sklearn.ensemble import VotingClassifier

from sklearn.model_selection import train_test_split
#from sklearn.metrics import accuracy_score,mean_absolute_error


#############

#list of descriptor label

coloumn_label_ls = ['group','period','H', 'C', 'N', 'O', 'B', 'F', 'Cl', 'Br', 'Si', 'P', 'S',
                 'ring_3','ring_4','ring_5','ring_6','ring_7','IsInRing', 
                  'CC_0', 'CC_1', 'CC_2', 'CC_3', 'CC_4', 'CC_5', 'CC_6', 'CC_7', 'CC_8', 
                  'CC_9', 'CN_0', 'CN_1', 'CN_2', 'CN_3', 'CN_4', 'CN_5', 'CN_6', 'CN_7', 
                  'CN_8', 'CN_9', 'CH_0', 'CH_1', 'CH_2', 'CH_3', 'CH_4', 'CH_5', 'CCl_0', 
                  'CCl_1', 'CCl_2', 'CCl_3', 'CCl_4', 'CCl_5', 'CO_0', 'CO_1', 'CO_2', 
                  'CO_3', 'CO_4', 'HN_0', 'HN_1', 'HN_2', 'HN_3', 'HN_4', 'HN_5', 'HN_6', 
                  'HN_7', 'HN_8', 'HN_9', 'CS_0', 'CS_1', 'CS_2', 'CS_3', 'CS_4', 'CS_5',
                  'CS_6', 'CS_7', 'CS_8', 'OS_0', 'OS_1', 'OS_2', 'OS_3', 'CF_0', 'CF_1', 
                  'CF_2', 'CF_3', 'CF_4', 'CF_5', 'CF_6', 'NO_0', 'NO_1', 'NO_2', 'NO_3', 
                  'NO_4', 'HO_0', 'HO_1', 'HO_2', 'HO_3', 'HO_4', 'HO_5', 'HO_6', 'HO_7',
                  'HO_8', 'HO_9', 'OP_0', 'OP_1', 'OP_2', 'OP_3',
                  'nH', 'nC', 'nN', 'nO', 'nB', 'nF', 'nCl', 'nBr', 'nSi', 'nP', 'nS', 'tn',
                  'nex1','n1H', 'n1C', 'n1N', 'n1O', 'n1B', 'n1F', 'n1Cl', 'n1Br', 'n1Si', 'n1P', 'n1S', 'tn1',
                  'nex2','n2H', 'n2C', 'n2N', 'n2O', 'n2B', 'n2F', 'n2Cl', 'n2Br', 'n2Si', 'n2P', 'n2S', 'tn2',
                  'nex3','n3H', 'n3C', 'n3N', 'n3O', 'n3B', 'n3F', 'n3Cl', 'n3Br', 'n3Si', 'n3P', 'n3S', 'tn3',
                  'nex4','n3H', 'n4C', 'n4N', 'n4O', 'n4B', 'n4F', 'n4Cl', 'n4Br', 'n4Si', 'n4P', 'n4S', 'tn4']


today = date.today()
date_str=today.strftime("%d%m%Y")


def get_LoopingLs(df, rxn_no='all', v=1):
    
    if rxn_no=='all':
        data_df=df.copy()
    else:
        data_df=df[df['code']<=rxn_no]
    
    
    unique_code_ls=data_df['code'].unique().tolist()
    
    
    corrs_idx_ls=[]
    for i in unique_code_ls:
        sub_df=data_df[data_df['code']==i]
        corrs_idx_ls.append(sub_df['idx'].tolist())
    
    n=len(unique_code_ls)

    full_ls=[i for i in range(0,n)]
    nested_list = []

    sub_idx_ls=list(range(0, len(full_ls)+1, v))

    for idx in range(0, len(sub_idx_ls)):
        try:
            nested_list.append(full_ls[sub_idx_ls[idx]:sub_idx_ls[idx+1]])
        except IndexError:
            print('error: idx='+str(full_ls[sub_idx_ls[idx]]))


    looping_ls=[]

    for idx_ls in nested_list:
        sub_ls=[i for i in full_ls if i not in idx_ls]
        looping_ls.append([sub_ls,idx_ls])

    return looping_ls


def list_duplicates(seq):
    ## from a list of items to [item, [list of indexes at which the item appears in the list]]
    ## e.g. input: [1, 2, 3, 4, 4, 3]
    ## output: [[1, [0]], [2, [1]], [3, [2, 5]], [4, [3, 4]]]

    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)  
    
    return [[key,locs] for key,locs in tally.items() if len(locs)>0]


#####
        


def to_at_ls(subls_ls):
    
    at_ls=[]
    for ls in subls_ls:
        for i in ls:
            at_ls.append(i)
            
    return np.array(at_ls)



class run_ML:
    def __init__(self,map_rxn_df):
        
        self.map_rxn_df = map_rxn_df
        
    
    def generate_data(self,method='m3',bl_method='rdkit',descriptor_arr ='2-bond+'):
        
        ## generating the descriptors 
        
        self.react_ls=[]

        self.r_descriptor_ls=[]
        
        self.r_label_ls=[]


        pass_fail_ls=[]
        error_ls=[]
        
        pass_sub_df_ls =[]
        
        code_ls=self.map_rxn_df['code'].unique().tolist()


        for code in code_ls:
            #print(idx)
            try:
        
                r_descriptor, r_label = tp.multi_path_descriptor_n_label(self.map_rxn_df, code, method=method, bl_method=bl_method,
                    descriptor_arr =descriptor_arr)
                
                
                self.r_label_ls.append(r_label)
                self.r_descriptor_ls.append(r_descriptor)
                
                
                pass_fail_ls.append('pass')
                error_ls.append('-')
                
                sub_df = self.map_rxn_df[self.map_rxn_df['code'] == code]
                pass_sub_df_ls.append(sub_df)
            
            except Exception as e:
                print('error: ',str(code))
                error_ls.append(str(e))
                pass_fail_ls.append('fail')
                print(e)
                
        
        self.status_df=pd.DataFrame({'idx':[i for i in range(0,len(code_ls))], 
                                'status': pass_fail_ls, 'error': error_ls})
        
        pass_no=len(self.status_df[self.status_df['status']=='pass'])
        self.pass_idx_ls=[code_ls[idx] for idx in range(0, len(pass_fail_ls)) if pass_fail_ls[idx]=='pass']
        idx_ls=[i for i in range(0,pass_no)]
        self.pass_idx_dict={i:j for i,j in zip(idx_ls, self.pass_idx_ls)}
        self.rev_pass_idx_dict={j:i for i,j in zip(idx_ls, self.pass_idx_ls)}
        self.pass_sub_df = pd.concat(pass_sub_df_ls)
        
        
        
        
                
                
    def split_data(self,by='rxn',incl='r',test_size=0.2, how='random'):
        
        ## split the data 
        ## how -- insert a list with this format [[training],[testing]] <-- should be the index not the actual descriptor or label
        
        self.by=by
        
        if incl=='r':
        
            descriptor_ls = self.r_descriptor_ls
            label_ls = self.r_label_ls
         
    
        #elif incl=='r+p':
        
            #descriptor_ls = self.r_descriptor_ls + self.p_descriptor_ls
            #label_ls = self.r_label_ls + self.p_label_ls
            
        

        
        if by=='atom':
            
            
            descriptor_arr=np.array(to_at_ls(descriptor_ls))
            label_arr=np.array(to_at_ls(label_ls))
            
            self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(descriptor_arr, label_arr, 
                                                            test_size=test_size, random_state=0)
            
            
        elif by=='rxn':
            
            idx_ls=[i for i in range(0,len(descriptor_ls))]
            
            if how=='random':

                self.idx_train, self.idx_test, idx_train, idx_test = train_test_split(idx_ls, idx_ls, test_size=test_size, random_state=50)
                
            else: 
                
                self.idx_train = how[0]
                self.idx_test = how[1]
                idx_train = how[0]
                idx_test = how[1]
                
            
            rxn_X_train = [descriptor_ls[idx] for idx in self.idx_train]
            self.rxn_X_test = [descriptor_ls[idx] for idx in self.idx_test]
            
            rxn_y_train = [label_ls[idx] for idx in self.idx_train]
            self.rxn_y_test = [label_ls[idx] for idx in self.idx_test]
            
            self.X_train = to_at_ls(rxn_X_train)
            self.y_train = to_at_ls(rxn_y_train)
            
    
    def training(self, model='RF', save_model = False, file_name='RF_model_',
        hyperpara={'n_estimators':100,'max_features':'sqrt','criterion':'gini'}):
        
        if model == 'RF':

            self.regr = RandomForestClassifier(n_estimators=hyperpara.get('n_estimators'),
                max_features=hyperpara.get('max_features'),criterion=hyperpara.get('criterion'))
            
        elif model == 'KNeigh':
            self.regr = KNeighborsClassifier(n_neighbors=5)
            
        elif model == 'SVC':
            self.regr = SVC()
        
        elif model == 'Gaussian':
            self.regr = GaussianProcessClassifier()
        
        elif model == 'NN':
            self.regr = MLPClassifier(solver='lbfgs', alpha=1e-4, random_state=None)
        
        
        if save_model == False:
        
            ## training the model 
         
            self.regr.fit(self.X_train, self.y_train)
            
        elif save_model == True:
            
            descriptor_ls = self.r_descriptor_ls
            label_ls = self.r_label_ls
            
            self.X_train=np.array(to_at_ls(descriptor_ls))
            self.y_train=np.array(to_at_ls(label_ls))
            
            
            self.regr.fit(self.X_train, self.y_train)
            
            
            filename = file_name+date_str+'.sav'
            pickle.dump(self.regr, open(filename, 'wb'))
            


    def evaluation(self, column=coloumn_label_ls):
        
        if self.by=='atom':
            
            y_predict=self.regr.predict(self.X_test)
            y_test=self.y_test

        elif self.by=='rxn':
            
            pred_label_ls=[]
            descriptor_ls=[]
            true_label_ls=[]
            idx_ls_ls=[]
            pred_label_ls_ls=[]
            
            atom_idx_ls=[]
            
            
            for idx in range(0,len(self.rxn_X_test)):
                sub_ls=[]

                for i in self.rxn_X_test[idx]:
                    at_y_predict=self.regr.predict([i])

                    sub_ls.append(at_y_predict[0])
                    descriptor_ls.append(i)
                    pred_label_ls_ls.append(at_y_predict[0])

                for j in self.rxn_y_test[idx]:
                    true_label_ls.append(j)
                    idx_ls_ls.append(idx)
                
                atom_idx_ls.extend([i for i in range(0,len(self.rxn_y_test[idx]))])
                
                
    
                pred_label_ls.append(sub_ls)



        rxn_idx_ls=[self.pass_idx_dict.get(i) for i in self.idx_test]
        long_rxn_idx_ls=[self.pass_idx_dict.get(i) for i in idx_ls_ls]
        
        self.label_df=pd.DataFrame({'idx': self.idx_test, 'rxn_idx':rxn_idx_ls, 
                                        'pred_label':pred_label_ls, 'actual_label':self.rxn_y_test}) 
        
        '''
        self.descriptor_df = pd.DataFrame({'rxn_idx':long_rxn_idx_ls,'atom_idx':atom_idx_ls ,
                                           'descriptor': descriptor_ls, 'true_label':true_label_ls, 
                                           'pred_label':pred_label_ls_ls})
    
    
        
        self.descriptor_label_df = pd.DataFrame(descriptor_ls, columns=coloumn_label_ls)
        
        '''
     
                
                
                
                
       
def get_label_ls(data_df):
    actual_label_ls=[]
    #print(data_df['actual_label'])

    for ls in data_df['actual_label'].tolist():
        actual_label_ls.extend(ls)
        
    pred_label_ls=[]

    for ls in data_df['pred_label'].tolist():
        pred_label_ls.extend(ls)
        
    return actual_label_ls,pred_label_ls



def get_metrice(pred_ls,test_ls):
    
    accuracy = accuracy_score(test_ls,pred_ls)
    precision = precision_score(test_ls,pred_ls)
    recall = recall_score(test_ls,pred_ls)

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
    
    True_label_pred=len([i for i in pred_ls if i >0])/len(pred_ls)
    True_label_act=len([i for i in test_ls if i >0])/len(test_ls)

    result_dict={'accuracy':accuracy,'precision':precision,'recall': recall,'tp':r_tp,'tn':r_tn, 
                'fp':r_fp, 'fn':r_fn,'postive_label_pred':True_label_pred,
                'postive_label_act':True_label_act}
    
    return result_dict
                
    


