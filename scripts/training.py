# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 07:05:52 2020

@author: hcji
"""

import pandas as pd
from chemvae.hyperparameters import load_params
from chemvae.train_vae import main_no_prop

from utils.utils import augment_smiles

compounds = pd.read_csv('data/ChEMBLv23_Dtraining_All_Scaffolds.txt', sep=' ')
dti = pd.read_csv('data/DEEPScreen_finalized_training_dataset_act_inact.txt', sep='\t', header=None)
dti.columns = ['target', 'drug']

params = load_params('params/params.json')

for i in dti.index:
    target = dti['target'][i]
    drugs = dti['drug'][i].split(',')
    drug_smiles = [compounds['SMILES'][list(compounds['Name']).index(d)] for d in drugs]
    aug_drug_smiles = augment_smiles(drug_smiles, num = 50)
    aug_drug_smiles = pd.DataFrame({'smiles':aug_drug_smiles})
    
    data_file = 'data/inputs/'+target+'.csv'
    aug_drug_smiles.to_csv(data_file, index=False)
    params['data_file'] = data_file
    main_no_prop(params)
    
    