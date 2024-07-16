import pandas as pd
from rdkit import Chem
from itertools import combinations

    
def get_inchi_ls(data1,with_product=False):
    
    rxns_ls = data1['reaction'].tolist()

    reactant=[i.split('>>')[0] for i in rxns_ls]

    if with_product == True:
        product=[i.split('>>')[1] for i in rxns_ls]

    inchi_rxn_ls=[]

    for idx in range(len(rxns_ls)):
        r_mol=Chem.MolFromSmiles(reactant[idx])
        r_inchi=Chem.MolToInchi(r_mol)

        if with_product == True:

            p_mol=Chem.MolFromSmiles(product[idx])
            p_inchi=Chem.MolToInchi(p_mol)

        inchi_rxn_ls.append(r_inchi+', '+p_inchi)
    
    return inchi_rxn_ls


file_ls=['./dataset/first-year_data_all_27122023.csv', './dataset/cyclo_data_v2_16072024.csv', './dataset/diels_alder_data_v7_19052024.csv',
         './dataset/fav_RGD1_13012024_v2.csv','./dataset/green_04012024.csv']


result = [list(comb) for comb in combinations(file_ls, 2)]


print('Test1: competitive reactions: reactions with the same reactants but different product found in both datasets')

for ls in result:
    print(ls)
    data1 = pd.read_csv(ls[0])
    data2 = pd.read_csv(ls[1])

    inchi_ls1=get_inchi_ls(data1)
    inchi_ls2=get_inchi_ls(data2)

    dup_ls=[]

    for idx in range(len(inchi_ls1)):
        if inchi_ls1[idx] in inchi_ls2:
            dup_ls.append(idx)

    print(dup_ls)

    if len(dup_ls) > 0:
        print('error')


print('Test2: duplicate reactions: the same reaction found in both datasets')

for ls in result:
    print(ls)
    data1 = pd.read_csv(ls[0])
    data2 = pd.read_csv(ls[1])

    inchi_ls1=get_inchi_ls(data1,with_product=True)
    inchi_ls2=get_inchi_ls(data2,with_product=True)

    dup_ls=[]

    for idx in range(len(inchi_ls1)):
        if inchi_ls1[idx] in inchi_ls2:
            dup_ls.append(idx)

    print(dup_ls)

    if len(dup_ls) > 0:
        print('error')


