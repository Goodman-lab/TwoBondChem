import pandas as pd 
from collections import defaultdict
from rdkit import Chem

def list_duplicates(seq):
    ## from a list of items to [item, [list of indexes at which the item appears in the list]]
    ## e.g. input: [1, 2, 3, 4, 4, 3]
    ## output: [[1, [0]], [2, [1]], [3, [2, 5]], [4, [3, 4]]]

    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)  
    
    return [[locs,key] for key,locs in tally.items() if len(locs)>0]


###########################
## RGD1 dataset initial processing 



cut_off=40

full_data=pd.read_csv('RGD1CHNO_AMsmiles.csv').rename(columns={"reaction": "Rind"})
fav_data=full_data[(full_data['DH']<0)&(full_data['DG_F']<cut_off)].copy()

react_ls=fav_data['reactant'].tolist()
prod_ls=fav_data['product'].tolist()
rxn_ls=[i+'>>'+j for i,j in zip(react_ls,prod_ls)]

inchi_react_ls=[]

for react in react_ls:
    r_mol=Chem.MolFromSmiles(react)
    r_inchi=Chem.MolToInchi(r_mol)
    inchi_react_ls.append(r_inchi)
sort_inchi_react_ls =list_duplicates(inchi_react_ls)

idx_matching_dict={}
for idx in range(0,len(sort_inchi_react_ls)):
    for i in sort_inchi_react_ls[idx][0]:
        idx_matching_dict[i]=idx
        
idx_ls=list(range(0,len(rxn_ls)))
code_ls=[idx_matching_dict.get(idx) for idx in idx_ls]

fav_data['code']=code_ls

fav_data['reaction']=rxn_ls

fav_data1=fav_data.sort_values(by='code')

idx_ls=[i for i in range(0,len(code_ls))]
fav_data2 = fav_data1[['code','reaction','DG_F','Rind']].copy()

fav_data2['idx']=idx_ls
fav_data3=fav_data2[['idx','code','reaction','DG_F','Rind']].copy()


fav_data3.to_csv('rgd_ea40.csv', index=False)

###########################
## Green dataset initial processing 


cut_off=40


full_data=pd.read_csv('b97d3.csv')
fav_data=full_data[(full_data['dh']<0)&(full_data['ea']<cut_off)].copy()


react_ls=fav_data['rsmi'].tolist()
prod_ls=fav_data['psmi'].tolist()
rxn_ls=[i+'>>'+j for i,j in zip(react_ls,prod_ls)]

inchi_react_ls=[]

for react in react_ls:
    r_mol=Chem.MolFromSmiles(react)
    r_inchi=Chem.MolToInchi(r_mol)
    inchi_react_ls.append(r_inchi)

sort_inchi_react_ls =list_duplicates(inchi_react_ls)

idx_matching_dict={}
for idx in range(0,len(sort_inchi_react_ls)):
    for i in sort_inchi_react_ls[idx][0]:
        idx_matching_dict[i]=idx
        
idx_ls=list(range(0,len(rxn_ls)))
code_ls=[idx_matching_dict.get(idx) for idx in idx_ls]

fav_data['code']=code_ls
fav_data['reaction']=rxn_ls

fav_data1=fav_data.sort_values(by='code')

idx_ls=[i for i in range(0,len(code_ls))]
fav_data2 = fav_data1[['code','reaction','ea']].copy()

fav_data2['idx']=idx_ls
fav_data3=fav_data2[['idx','code','reaction','ea']].copy()

fav_data3.to_csv('green_ea40.csv',  index=False)


###########################
## Diels-Alder dataset initial processing 

full_data=pd.read_csv('diels_alder_data.csv')


rxn_ls1=full_data['reaction'].tolist()
rxn_ls=[]
for i in rxn_ls1:
    if '*' not in i:
        rxn_ls.append(i)
        
        
react_ls=[]
prod_ls=[]
for rxn in rxn_ls:
    sub_ls=rxn.split('>>')
    react_ls.append(sub_ls[0])
    prod_ls.append(sub_ls[1])

inchi_react_ls=[]

for react in react_ls:
    r_mol=Chem.MolFromSmiles(react)
    r_inchi=Chem.MolToInchi(r_mol)
    inchi_react_ls.append(r_inchi)

sort_inchi_react_ls =list_duplicates(inchi_react_ls)

idx_matching_dict={}
for idx in range(0,len(sort_inchi_react_ls)):
    for i in sort_inchi_react_ls[idx][0]:
        idx_matching_dict[i]=idx
        
idx_ls=list(range(0,len(rxn_ls)))
code_ls=[idx_matching_dict.get(idx) for idx in idx_ls]


new_data_df=pd.DataFrame({'reaction':rxn_ls})
new_data_df['code']=code_ls

full_data2=new_data_df.sort_values('code').copy()


full_data2['idx']=idx_ls


full_data3=full_data2[['idx','code','reaction']]

full_data3.to_csv('diels_alder_data_v3.csv',index=False)


###########################
## [3+2] cycloaddition dataset initial processing 


full_data=pd.read_csv('reactions_combined_finalized.csv')


rxn_ls=full_data['reaction'].tolist()
        
        
react_ls=[]
prod_ls=[]
for rxn in rxn_ls:
    sub_ls=rxn.split('>>')
    react_ls.append(sub_ls[0])
    prod_ls.append(sub_ls[1])

inchi_react_ls=[]

for react in react_ls:
    r_mol=Chem.MolFromSmiles(react)
    r_inchi=Chem.MolToInchi(r_mol)
    inchi_react_ls.append(r_inchi)

sort_inchi_react_ls =list_duplicates(inchi_react_ls)

idx_matching_dict={}
for idx in range(0,len(sort_inchi_react_ls)):
    for i in sort_inchi_react_ls[idx][0]:
        idx_matching_dict[i]=idx
        
idx_ls=list(range(0,len(rxn_ls)))
code_ls=[idx_matching_dict.get(idx) for idx in idx_ls]


new_data_df=pd.DataFrame({'reaction':rxn_ls})
new_data_df['code']=code_ls

full_data2=new_data_df.sort_values('code').copy()


full_data2['idx']=idx_ls


full_data3=full_data2[['idx','code','reaction']]


full_data3.to_csv('cyclo.csv',index=False)


