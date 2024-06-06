#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 13:30:18 2023

@author: chingchinglam
"""

from datetime import date
import os

#from rdkit import Chem
from atidx import AtIdx,mapping_competitive_pathways_all
import pandas as pd


file_name = 'exam_test.csv'
reaction_df = pd.read_csv(file_name)
result=mapping_competitive_pathways_all(reaction_df)

today = date.today()
date_str=today.strftime("%d%m%Y")
result.to_csv(file_name[:-4]+'_'+date_str+'.csv' )


