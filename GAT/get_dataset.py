from typing import List, Tuple, Union
from itertools import zip_longest
import logging
from torch_geometric.data import Data, DataLoader, InMemoryDataset
from rdkit import Chem
import torch
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from chemprop.rdkit import make_mol
from tqdm import tqdm
import torch
from torch.utils.data import Dataset

class Featurization_parameters:
    """
    A class holding molecule featurization parameters as attributes.
    """
    def __init__(self) -> None:

        # Atom feature sizes
        self.MAX_ATOMIC_NUM = 100
        self.ATOM_FEATURES = {
            'atomic_num': list(range(self.MAX_ATOMIC_NUM)),
            'degree': [0, 1, 2, 3, 4, 5],
            'formal_charge': [-1, -2, 1, 2, 0],
            'chiral_tag': [0, 1, 2, 3],
            'num_Hs': [0, 1, 2, 3, 4],
            'hybridization': [
                Chem.rdchem.HybridizationType.SP,
                Chem.rdchem.HybridizationType.SP2,
                Chem.rdchem.HybridizationType.SP3,
                Chem.rdchem.HybridizationType.SP3D,
                Chem.rdchem.HybridizationType.SP3D2
            ],
        }

        # Distance feature sizes
        self.PATH_DISTANCE_BINS = list(range(10))
        self.THREE_D_DISTANCE_MAX = 20
        self.THREE_D_DISTANCE_STEP = 1
        self.THREE_D_DISTANCE_BINS = list(range(0, self.THREE_D_DISTANCE_MAX + 1, self.THREE_D_DISTANCE_STEP))

        # len(choices) + 1 to include room for uncommon values; + 2 at end for IsAromatic and mass
        self.ATOM_FDIM = sum(len(choices) + 1 for choices in self.ATOM_FEATURES.values()) + 2
        self.EXTRA_ATOM_FDIM = 0
        self.BOND_FDIM = 14
        self.EXTRA_BOND_FDIM = 0
        self.REACTION_MODE = None
        self.EXPLICIT_H = False
        self.REACTION = False
        self.ADDING_H = False
        self.KEEP_ATOM_MAP = False

# Create a global parameter object for reference throughout this module
PARAMS = Featurization_parameters()

def bond_features(bond: Chem.rdchem.Bond) -> List[Union[bool, int, float]]:
    """
    Builds a feature vector for a bond.

    :param bond: An RDKit bond.
    :return: A list containing the bond features.
    """
    if bond is None:
        fbond = [1] + [0] * (PARAMS.BOND_FDIM - 1)
    else:
        bt = bond.GetBondType()
        fbond = [
            0,  # bond is not None
            bt == Chem.rdchem.BondType.SINGLE,
            bt == Chem.rdchem.BondType.DOUBLE,
            bt == Chem.rdchem.BondType.TRIPLE,
            bt == Chem.rdchem.BondType.AROMATIC,
            (bond.GetIsConjugated() if bt is not None else 0),
            (bond.IsInRing() if bt is not None else 0)
        ]
        fbond += onek_encoding_unk(int(bond.GetStereo()), list(range(6)))
    return fbond

def onek_encoding_unk(value: int, choices: List[int]) -> List[int]:
    """
    Creates a one-hot encoding with an extra category for uncommon values.

    :param value: The value for which the encoding should be one.
    :param choices: A list of possible values.
    :return: A one-hot encoding of the :code:`value` in a list of length :code:`len(choices) + 1`.
             If :code:`value` is not in :code:`choices`, then the final element in the encoding is 1.
    """
    encoding = [0] * (len(choices) + 1)
    index = choices.index(value) if value in choices else -1
    encoding[index] = 1

    return encoding


def atom_features(atom: Chem.rdchem.Atom, functional_groups: List[int] = None) -> List[Union[bool, int, float]]:
    """
    Builds a feature vector for an atom.

    :param atom: An RDKit atom.
    :param functional_groups: A k-hot vector indicating the functional groups the atom belongs to.
    :return: A list containing the atom features.
    """
    if atom is None:
        features = [0] * PARAMS.ATOM_FDIM
    else:
        features = onek_encoding_unk(atom.GetAtomicNum() - 1, PARAMS.ATOM_FEATURES['atomic_num']) + \
            onek_encoding_unk(atom.GetTotalDegree(), PARAMS.ATOM_FEATURES['degree']) + \
            onek_encoding_unk(atom.GetFormalCharge(), PARAMS.ATOM_FEATURES['formal_charge']) + \
            onek_encoding_unk(int(atom.GetChiralTag()), PARAMS.ATOM_FEATURES['chiral_tag']) + \
            onek_encoding_unk(int(atom.GetTotalNumHs()), PARAMS.ATOM_FEATURES['num_Hs']) + \
            onek_encoding_unk(int(atom.GetHybridization()), PARAMS.ATOM_FEATURES['hybridization']) + \
            [1 if atom.GetIsAromatic() else 0] + \
            [atom.GetMass() * 0.01]  # scaled to about the same range as other features
        if functional_groups is not None:
            features += functional_groups
    return features

def atom_features_zeros(atom: Chem.rdchem.Atom) -> List[Union[bool, int, float]]:
    """
    Builds a feature vector for an atom containing only the atom number information.

    :param atom: An RDKit atom.
    :return: A list containing the atom features.
    """
    if atom is None:
        features = [0] * PARAMS.ATOM_FDIM
    else:
        features = onek_encoding_unk(atom.GetAtomicNum() - 1, PARAMS.ATOM_FEATURES['atomic_num']) + \
            [0] * (PARAMS.ATOM_FDIM - PARAMS.MAX_ATOMIC_NUM - 1) #set other features to zero
    return features

def map_reac_to_prod(mol_reac: Chem.Mol, mol_prod: Chem.Mol):
    """
    Build a dictionary of mapping atom indices in the reactants to the products.

    :param mol_reac: An RDKit molecule of the reactants.
    :param mol_prod: An RDKit molecule of the products.
    :return: A dictionary of corresponding reactant and product atom indices.
    """
    only_prod_ids = []
    prod_map_to_id = {}
    mapnos_reac = set([atom.GetAtomMapNum() for atom in mol_reac.GetAtoms()]) 
    for atom in mol_prod.GetAtoms():
        mapno = atom.GetAtomMapNum()
        if mapno > 0:
            prod_map_to_id[mapno] = atom.GetIdx()
            if mapno not in mapnos_reac:
                only_prod_ids.append(atom.GetIdx()) 
        else:
            only_prod_ids.append(atom.GetIdx())
    only_reac_ids = []
    reac_id_to_prod_id = {}
    for atom in mol_reac.GetAtoms():
        mapno = atom.GetAtomMapNum()
        if mapno > 0:
            try:
                reac_id_to_prod_id[atom.GetIdx()] = prod_map_to_id[mapno]
            except KeyError:
                only_reac_ids.append(atom.GetIdx())
        else:
            only_reac_ids.append(atom.GetIdx())
    return reac_id_to_prod_id, only_prod_ids, only_reac_ids


def get_rxn_info_from_sma(smarts):
    mol = smarts.split('>>')
    mol_reac = Chem.MolFromSmiles(mol[0])
    mol_prod = Chem.MolFromSmiles(mol[1])

    # Atom features
    ri2pi, pio, rio = map_reac_to_prod(mol_reac, mol_prod)
    f_atoms_reac = [atom_features(atom) for atom in mol_reac.GetAtoms()] + [atom_features_zeros(mol_prod.GetAtomWithIdx(index)) for index in pio]
    f_atoms_prod = [atom_features(mol_prod.GetAtomWithIdx(ri2pi[atom.GetIdx()])) if atom.GetIdx() not in rio else
                                atom_features_zeros(atom) for atom in mol_reac.GetAtoms()] + [atom_features(mol_prod.GetAtomWithIdx(index)) for index in pio]
    f_atoms_diff = [list(map(lambda x, y: x - y, ii, jj)) for ii, jj in zip(f_atoms_prod, f_atoms_reac)]
    f_atoms = [x+y[PARAMS.MAX_ATOMIC_NUM+1:] for x,y in zip(f_atoms_reac, f_atoms_diff)]
    n_atoms = len(f_atoms)
    n_atoms_reac = mol_reac.GetNumAtoms()

    # Bond features
    f_bonds = []
    edge_dir1 = []
    edge_dir2 = []
    for a1 in range(n_atoms):
        for a2 in range(a1 + 1, n_atoms):
            if a1 >= n_atoms_reac and a2 >= n_atoms_reac: # Both atoms only in 
                bond_prod = mol_prod.GetBondBetweenAtoms(pio[a1 - n_atoms_reac], pio[a2 - n_atoms_reac])
                bond_reac = None
            elif a1 < n_atoms_reac and a2 >= n_atoms_reac: # One atom only in product
                bond_reac = None
                if a1 in ri2pi.keys():
                    bond_prod = mol_prod.GetBondBetweenAtoms(ri2pi[a1], pio[a2 - n_atoms_reac])
                else:
                    bond_prod = None # Atom atom only in reactant, the other only in 
            else:
                bond_reac = mol_reac.GetBondBetweenAtoms(a1, a2)
                if a1 in ri2pi.keys() and a2 in ri2pi.keys():
                    bond_prod = mol_prod.GetBondBetweenAtoms(ri2pi[a1], ri2pi[a2]) #Both atoms in both reactant and product
                else:  
                    bond_prod = None # One or both atoms only in reactant
            if bond_reac is None and bond_prod is None:
                continue

            f_bond_reac = bond_features(bond_reac)
            f_bond_prod = bond_features(bond_prod)
            f_bond_diff = [y - x for x, y in zip(f_bond_reac, f_bond_prod)]
            
            f_bond = f_bond_reac + f_bond_diff
            f_bonds.append(f_atoms[a1] + f_bond)
            f_bonds.append(f_atoms[a2] + f_bond)
            edge_dir1.append(a1)
            edge_dir2.append(a2)
            edge_dir1.append(a2)
            edge_dir2.append(a1)
            edge_2d = torch.from_numpy(np.array([edge_dir1, edge_dir2]))
    n_atoms = torch.tensor([n_atoms], dtype=torch.long)
    f_atoms = torch.tensor(f_atoms, dtype=torch.float)
    f_bonds = torch.tensor(f_bonds, dtype=torch.float)
    return n_atoms,f_atoms, f_bonds, edge_2d

def process_single_graph(smarts, target, index):
    n_atoms, f_atoms, f_bonds, edge_index = get_rxn_info_from_sma(smarts)
    batch = torch.full((n_atoms,), index, dtype=torch.long)
    if target != None:
        y = torch.tensor([target], dtype=torch.long)
    else:
        y = y = torch.tensor([0], dtype=torch.long)
    return n_atoms, f_atoms, f_bonds, edge_index, y, batch


class GraphDataset(Dataset):
    def __init__(self, data, target):
        self.target = target
        self.data = data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        smart = self.data['reaction'][idx]
        if self.target != None:
            target = self.data[self.target][idx]
        else:
            target = None
        n_atoms, f_atoms, f_bonds, edge_index, y, batch = process_single_graph(smart, target, idx)
        
        return {
            'x': f_atoms,
            'edge_index': edge_index,
            'edge_attr': f_bonds,
            'y': y,
            'batch': batch
        }
    

class GraphDataLoader:
    def __init__(self, dataset , batch_size: int, collate_fn):
        self.dataset = dataset
        self.batch_size = batch_size
        self.collate_fn = collate_fn

    def __iter__(self):
        self.current_index = 0
        return self

    def __next__(self):
        if self.current_index >= len(self.dataset):
            raise StopIteration
        
        batch = []
        for i in range(self.batch_size):
            if self.current_index + i >= len(self.dataset):
                break
            batch.append(self.dataset[self.current_index + i])
        self.current_index += self.batch_size
        
        return self.collate_fn(batch)
    
    def __len__(self):
        return len(self.dataset) // self.batch_size + 1

def collate_fn(data):
    data_x = [d['x'] for d in data]  # Node features
    data_edge_index = [d['edge_index'] for d in data]  # Edge indices
    data_edge_attr = [d['edge_attr'] for d in data]  # Edge attributes
    data_y = [d['y'] for d in data]  # Targets
    # Assuming data_batch was to mark which graph a node belongs to, we'll recalculate it here
    # data_batch = [d['batch'] for d in data]

    # Concatenate node features
    data_x = torch.cat(data_x, dim=0)
    # Concatenate edge attributes
    data_edge_attr = torch.cat(data_edge_attr, dim=0)
    # Concatenate targets
    data_y = torch.cat(data_y, dim=0)

    # Handle edge_index
    cumsum_node = 0  # Cumulative number of nodes
    data_edge_index_new = []
    for edge_index in data_edge_index:
        # Update edge_index with the cumulative sum to adjust indices for the merged graph
        edge_index_updated = edge_index + cumsum_node
        data_edge_index_new.append(edge_index_updated)
        # Update the cumulative node count
        cumsum_node += edge_index.max() + 1  # Assuming edge_index starts from 0 for each graph

    # Concatenate all updated edge_index tensors
    data_edge_index = torch.cat(data_edge_index_new, dim=1)

    # Handle batch information, i.e., to which graph each node belongs
    data_batch = []
    for i, d in enumerate(data):
        batch_size = d['x'].size(0)
        data_batch.append(torch.full((batch_size,), i, dtype=torch.long))
    data_batch = torch.cat(data_batch, dim=0)


    return Data(x=data_x, edge_index=data_edge_index, edge_attr=data_edge_attr, y=data_y, batch=data_batch)


if __name__ == '__main__':  
    pass
