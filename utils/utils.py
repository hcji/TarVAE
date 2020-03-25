# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 09:40:13 2020

@author: hcji
"""

import random
import rdkit.Chem as rkc
from tqdm import tqdm


def randomize_smiles(mol, random_type="restricted"):
    """
    Returns a random SMILES given a SMILES of a molecule.
    :param mol: A Mol object
    :param random_type: The type (unrestricted, restricted) of randomization performed.
    :return : A random SMILES string of the same molecule or None if the molecule is invalid.
    """
    if not mol:
        return None

    if random_type == "unrestricted":
        return rkc.MolToSmiles(mol, canonical=False, doRandom=True, isomericSmiles=False)
    if random_type == "restricted":
        new_atom_order = list(range(mol.GetNumAtoms()))
        random.shuffle(new_atom_order)
        random_mol = rkc.RenumberAtoms(mol, newOrder=new_atom_order)
        return rkc.MolToSmiles(random_mol, canonical=False, isomericSmiles=False)
    raise ValueError("Type '{}' is not valid".format(random_type))
    
    
def augment_smiles(smiles_list, num=10, random_type="restricted"):
    output = smiles_list.copy()
    for smi in tqdm(smiles_list):
        mol = rkc.MolFromSmiles(smi)
        if not mol:
            continue
        aug_smiles = [randomize_smiles(mol, random_type=random_type)for i in range(1, num)]
        output += aug_smiles
    return output


    
    