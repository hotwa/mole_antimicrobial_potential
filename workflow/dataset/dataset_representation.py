import os
import time
import yaml
import numpy as np
import pandas as pd
from contextlib import contextmanager
import threading

import torch
from torch_geometric.data import Data, Dataset, Batch
from torch_geometric.loader import DataLoader

from rdkit import Chem
from rdkit.Chem.rdchem import HybridizationType
from rdkit.Chem.rdchem import BondType as BT
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit import RDLogger                           
from rdkit.Chem.SaltRemover import SaltRemover    

from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import MolSurf
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors3D

RDLogger.DisableLog('rdApp.*')  


ATOM_LIST = list(range(1,119))
CHIRALITY_LIST = [
    Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
    Chem.rdchem.ChiralType.CHI_OTHER
]
BOND_LIST = [
    BT.SINGLE, 
    BT.DOUBLE, 
    BT.TRIPLE, 
    BT.AROMATIC
]
BONDDIR_LIST = [
    Chem.rdchem.BondDir.NONE,
    Chem.rdchem.BondDir.ENDUPRIGHT,
    Chem.rdchem.BondDir.ENDDOWNRIGHT
]
ATOM_INDEX = {atomic_num: index for index, atomic_num in enumerate(ATOM_LIST)}
CHIRALITY_INDEX = {tag: index for index, tag in enumerate(CHIRALITY_LIST)}
BOND_INDEX = {bond_type: index for index, bond_type in enumerate(BOND_LIST)}
BONDDIR_INDEX = {bond_dir: index for index, bond_dir in enumerate(BONDDIR_LIST)}
_DETERMINISTIC_TORCH_LOCK = threading.Lock()

# A FUNCTION TO READ SMILES from file 
def read_smiles(data_path, smile_col="rdkit_no_salt", id_col="prestwick_ID"):

    """
    Read SMILES data from a file and remove invalid SMILES.

    Parameters:
    - data_path (str): Path to the file containing SMILES data.
    - smile_col (str, optional): Name of the column containing SMILES strings (default is "rdkit_no_salt").
    - id_col (str, optional): Name of the column containing molecule IDs (default is "prestwick_ID").

    Returns:
    - smile_df (pandas.DataFrame): DataFrame containing SMILES data with specified columns.
    """
    
    # Read the data
    smile_df = pd.read_csv(data_path, sep='\t')
    smile_df = smile_df[[smile_col, id_col]]

    # Remove NaN
    smile_df = smile_df.dropna()

    # Remove invalid smiles
    smile_df = smile_df[smile_df[smile_col].apply(lambda x: Chem.MolFromSmiles(x) is not None)]

    return smile_df

# A FUNCTION TO CALCULATE DESCRIPTORS
def calc_descriptors(smiles, id):
    """
    Calculate molecular descriptors for a given SMILES string.
    ORIGINAL FUNCTION FROM ALGAVI & BORENSTEIN, 2023

    Parameters:
    - smiles (str): SMILES string of the molecule.
    - id (str): ID of the molecule.

    Returns:
    - smiles_df (pandas.DataFrame): DataFrame containing calculated molecular descriptors.
    """

    #define mol
    mol = Chem.MolFromSmiles(smiles)
    mol1=Chem.AddHs(mol)
    
    #remove salts
    remover = SaltRemover()
    mol1 =  remover.StripMol(mol1)
    smiles_df = pd.DataFrame({"chem_id": [id]})
    
    #calculate conformers
    AllChem.EmbedMolecule(mol1)
    AllChem.MMFFOptimizeMolecule(mol1)
    
    ## genral stracture descriptors (rdkit.Chem.Descriptors module)
    smiles_df["MolWt"] = Chem.Descriptors.MolWt(mol1)
    smiles_df["BertzCT"] = Chem.Descriptors.BertzCT(mol1)
    smiles_df["MolLogP"] = Chem.Descriptors.MolLogP(mol1)
    smiles_df["MolMR"] = Chem.Descriptors.MolMR(mol1)
    smiles_df["HeavyAtomCount"] = Chem.Descriptors.HeavyAtomCount(mol1)
    smiles_df["NumHAcceptors"] = Chem.Descriptors.NumHAcceptors(mol1)
    smiles_df["NumHDonors"] = Chem.Descriptors.NumHDonors(mol1)
    smiles_df["NumValenceElectrons"] = Chem.Descriptors.NumValenceElectrons(mol1)
    smiles_df["RingCount"] = Chem.Descriptors.RingCount(mol1)
    smiles_df["FractionCSP3"] = Chem.Descriptors.FractionCSP3(mol1)
    smiles_df["NHOHCount"] = Chem.Descriptors.NHOHCount(mol1)
    smiles_df["NOCount"] = Chem.Descriptors.NOCount(mol1)
    smiles_df["HeavyAtomMolWt"] = Chem.Descriptors.HeavyAtomMolWt(mol1)
    smiles_df["MaxAbsPartialCharge"] = Chem.Descriptors.MaxAbsPartialCharge(mol1)
    smiles_df["MaxPartialCharge"] = Chem.Descriptors.MaxPartialCharge(mol1)
    smiles_df["MinAbsPartialCharge"] = Chem.Descriptors.MinAbsPartialCharge(mol1)
    smiles_df["MinPartialCharge"] = Chem.Descriptors.MinPartialCharge(mol1)

    #Graph descriptors from Chem.rdMolDescriptors module
    smiles_df["Chi0n"] =Chem.rdMolDescriptors.CalcChi0n(mol1)
    smiles_df["Chi0v"] =Chem.rdMolDescriptors.CalcChi0v(mol1)
    smiles_df["Chi1n"] =Chem.rdMolDescriptors.CalcChi1n(mol1)
    smiles_df["Chi1v"] =Chem.rdMolDescriptors.CalcChi1v(mol1)
    smiles_df["Chi2n"] =Chem.rdMolDescriptors.CalcChi2n(mol1)
    smiles_df["Chi2v"] =Chem.rdMolDescriptors.CalcChi2v(mol1)
    smiles_df["Chi3n"] =Chem.rdMolDescriptors.CalcChi3n(mol1)
    smiles_df["Chi3v"] =Chem.rdMolDescriptors.CalcChi3v(mol1)
    smiles_df["Chi4n"] =Chem.rdMolDescriptors.CalcChi4n(mol1)
    smiles_df["Chi4v"] =Chem.rdMolDescriptors.CalcChi4v(mol1)
    smiles_df["HallKierAlpha"] =Chem.rdMolDescriptors.CalcHallKierAlpha(mol1)
    smiles_df["Kappa1"] =Chem.rdMolDescriptors.CalcKappa1(mol1)
    smiles_df["Kappa2"] =Chem.rdMolDescriptors.CalcKappa2(mol1)
    smiles_df["Kappa3"] =Chem.rdMolDescriptors.CalcKappa3(mol1)
    smiles_df["LabuteASA"] =Chem.rdMolDescriptors.CalcLabuteASA(mol1)
    smiles_df["NumAliphaticCarbocycles"] =Chem.rdMolDescriptors.CalcNumAliphaticCarbocycles(mol1)
    smiles_df["NumAliphaticHeterocycles"] =Chem.rdMolDescriptors.CalcNumAliphaticHeterocycles(mol1)
    smiles_df["NumAliphaticRings"] =Chem.rdMolDescriptors.CalcNumAliphaticRings(mol1)
    smiles_df["NumAmideBonds"] =Chem.rdMolDescriptors.CalcNumAmideBonds(mol1)
    smiles_df["NumAromaticCarbocycles"] =Chem.rdMolDescriptors.CalcNumAromaticCarbocycles(mol1)
    smiles_df["NumAromaticHeterocycles"] =Chem.rdMolDescriptors.CalcNumAromaticHeterocycles(mol1)
    smiles_df["NumAromaticRings"] =Chem.rdMolDescriptors.CalcNumAromaticRings(mol1)
    smiles_df["NumBridgeheadAtoms"] =Chem.rdMolDescriptors.CalcNumBridgeheadAtoms(mol1)
    smiles_df["NumHBA"] =Chem.rdMolDescriptors.CalcNumHBA(mol1)
    smiles_df["NumHBD"] =Chem.rdMolDescriptors.CalcNumHBD(mol1)
    smiles_df["NumHeteroatoms"] =Chem.rdMolDescriptors.CalcNumHeteroatoms(mol1)
    smiles_df["NumHeterocycles"] =Chem.rdMolDescriptors.CalcNumHeterocycles(mol1)
    smiles_df["NumLipinskiHBA"] =Chem.rdMolDescriptors.CalcNumLipinskiHBA(mol1)
    smiles_df["NumLipinskiHBD"] =Chem.rdMolDescriptors.CalcNumLipinskiHBD(mol1)
    smiles_df["NumRings"] =Chem.rdMolDescriptors.CalcNumRings(mol1)
    smiles_df["NumRotatableBonds"] =Chem.rdMolDescriptors.CalcNumRotatableBonds(mol1)
    smiles_df["NumSaturatedCarbocycles"] =Chem.rdMolDescriptors.CalcNumSaturatedCarbocycles(mol1)
    smiles_df["NumSaturatedHeterocycles"] =Chem.rdMolDescriptors.CalcNumSaturatedHeterocycles(mol1)
    smiles_df["NumSaturatedRings"] =Chem.rdMolDescriptors.CalcNumSaturatedRings(mol1)
    smiles_df["NumSpiroAtoms"] =Chem.rdMolDescriptors.CalcNumSpiroAtoms(mol1)
    smiles_df["PBF"] =Chem.rdMolDescriptors.CalcPBF(mol1)
    smiles_df["PMI1"] =Chem.rdMolDescriptors.CalcPMI1(mol1)
    smiles_df["PMI2"] =Chem.rdMolDescriptors.CalcPMI2(mol1)
    smiles_df["PMI3"] =Chem.rdMolDescriptors.CalcPMI3(mol1)
    smiles_df["TPSA"] =Chem.rdMolDescriptors.CalcTPSA(mol1)

    #surface properties from  
    ##PEOE =  The Partial Equalization of Orbital Electronegativities method of calculating atomic partial charges
    smiles_df["PEOE_VAS1"]=Chem.MolSurf.PEOE_VSA1(mol1)
    smiles_df["PEOE_VAS2"]=Chem.MolSurf.PEOE_VSA2(mol1)
    smiles_df["PEOE_VAS3"]=Chem.MolSurf.PEOE_VSA3(mol1)
    smiles_df["PEOE_VAS4"]=Chem.MolSurf.PEOE_VSA4(mol1)
    smiles_df["PEOE_VAS5"]=Chem.MolSurf.PEOE_VSA5(mol1)
    smiles_df["PEOE_VAS6"]=Chem.MolSurf.PEOE_VSA6(mol1)
    smiles_df["PEOE_VAS7"]=Chem.MolSurf.PEOE_VSA7(mol1)
    smiles_df["PEOE_VAS8"]=Chem.MolSurf.PEOE_VSA8(mol1)
    smiles_df["PEOE_VAS9"]=Chem.MolSurf.PEOE_VSA9(mol1)
    smiles_df["PEOE_VAS10"]=Chem.MolSurf.PEOE_VSA10(mol1)
    smiles_df["PEOE_VAS11"]=Chem.MolSurf.PEOE_VSA11(mol1)
    smiles_df["PEOE_VAS12"]=Chem.MolSurf.PEOE_VSA12(mol1)
    smiles_df["PEOE_VAS13"]=Chem.MolSurf.PEOE_VSA13(mol1)
    smiles_df["PEOE_VAS14"]=Chem.MolSurf.PEOE_VSA14(mol1)
    ##SMR = Molecular refractivity
    smiles_df["SMR_VSA1"]=Chem.MolSurf.SMR_VSA1(mol1)
    smiles_df["SMR_VSA2"]=Chem.MolSurf.SMR_VSA2(mol1)
    smiles_df["SMR_VSA3"]=Chem.MolSurf.SMR_VSA3(mol1)
    smiles_df["SMR_VSA4"]=Chem.MolSurf.SMR_VSA4(mol1)
    smiles_df["SMR_VSA5"]=Chem.MolSurf.SMR_VSA5(mol1)
    smiles_df["SMR_VSA6"]=Chem.MolSurf.SMR_VSA6(mol1)
    smiles_df["SMR_VSA7"]=Chem.MolSurf.SMR_VSA7(mol1)
    smiles_df["SMR_VSA8"]=Chem.MolSurf.SMR_VSA8(mol1)
    smiles_df["SMR_VSA9"]=Chem.MolSurf.SMR_VSA9(mol1)
    smiles_df["SMR_VSA10"]=Chem.MolSurf.SMR_VSA10(mol1)
    ##slogp = Log of the octanol/water partition coefficient
    smiles_df["SlogP_VSA1"]=Chem.MolSurf.SlogP_VSA1(mol1)
    smiles_df["SlogP_VSA2"]=Chem.MolSurf.SlogP_VSA2(mol1)
    smiles_df["SlogP_VSA3"]=Chem.MolSurf.SlogP_VSA3(mol1)
    smiles_df["SlogP_VSA4"]=Chem.MolSurf.SlogP_VSA4(mol1)
    smiles_df["SlogP_VSA5"]=Chem.MolSurf.SlogP_VSA5(mol1)
    smiles_df["SlogP_VSA6"]=Chem.MolSurf.SlogP_VSA6(mol1)
    smiles_df["SlogP_VSA7"]=Chem.MolSurf.SlogP_VSA7(mol1)
    smiles_df["SlogP_VSA8"]=Chem.MolSurf.SlogP_VSA8(mol1)
    smiles_df["SlogP_VSA9"]=Chem.MolSurf.SlogP_VSA9(mol1)
    smiles_df["SlogP_VSA10"]=Chem.MolSurf.SlogP_VSA10(mol1)
    smiles_df["SlogP_VSA11"]=Chem.MolSurf.SlogP_VSA11(mol1)
    ##others
    smiles_df["pyLabuteASA"]=Chem.MolSurf.pyLabuteASA(mol1)

    #3D descriptors from Chem.Descriptors3D module
    smiles_df["Asphericity"] = Chem.Descriptors3D.Asphericity(mol1)
    smiles_df["Eccentricity"] = Chem.Descriptors3D.Eccentricity(mol1)
    smiles_df["InertialShapeFactor"] = Chem.Descriptors3D.InertialShapeFactor(mol1)
    smiles_df["RadiusOfGyration"] = Chem.Descriptors3D.RadiusOfGyration(mol1)
    smiles_df["SpherocityIndex"] = Chem.Descriptors3D.SpherocityIndex(mol1)
    
    return(smiles_df)

# Here we can add more molecular descriptors
class MoleculeDataset(Dataset):

    """
    Dataset class for creating molecular graphs.

    Attributes:
    - smile_df (pandas.DataFrame): DataFrame containing SMILES data.
    - smile_column (str): Name of the column containing SMILES strings.
    - id_column (str): Name of the column containing molecule IDs.
    """

    def __init__(self, smile_df, smile_column, id_column, enable_profiling=False):
        super(Dataset, self).__init__()

        # Gather the SMILES and the corresponding IDs
        self.smiles_data = smile_df[smile_column].tolist()
        self.id_data = smile_df[id_column].tolist()
        self.enable_profiling = enable_profiling
        self._graph_cache = {}

    def _build_graph_payload(self, smiles):
        parse_start = time.perf_counter() if self.enable_profiling else None
        mol = Chem.MolFromSmiles(smiles)
        parse_end = time.perf_counter() if self.enable_profiling else None

        add_hs_start = time.perf_counter() if self.enable_profiling else None
        mol = Chem.AddHs(mol)
        add_hs_end = time.perf_counter() if self.enable_profiling else None

        atom_count = mol.GetNumAtoms()
        atom_features_np = np.empty((atom_count, 2), dtype=np.int64)
        atom_type_column = atom_features_np[:, 0]
        atom_chirality_column = atom_features_np[:, 1]
        has_unknown_atom = False

        atom_feature_start = time.perf_counter() if self.enable_profiling else None
        for atom_index, atom in enumerate(mol.GetAtoms()):
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 0:
                has_unknown_atom = True

            atom_type_column[atom_index] = ATOM_INDEX[atomic_num]
            atom_chirality_column[atom_index] = CHIRALITY_INDEX[atom.GetChiralTag()]
        atom_feature_end = time.perf_counter() if self.enable_profiling else None
        x = torch.from_numpy(atom_features_np)

        bond_count = mol.GetNumBonds()
        edge_index = torch.empty((2, bond_count * 2), dtype=torch.long)
        if bond_count == 0:
            edge_attr = torch.empty((0,), dtype=torch.long)
        else:
            edge_index_np = np.empty((2, bond_count * 2), dtype=np.int64)
            edge_attr_np = np.empty((bond_count * 2, 2), dtype=np.int64)
            edge_index0 = edge_index_np[0]
            edge_index1 = edge_index_np[1]
            edge_attr0 = edge_attr_np[:, 0]
            edge_attr1 = edge_attr_np[:, 1]

        bond_feature_start = time.perf_counter() if self.enable_profiling else None
        edge_offset = 0
        for bond in mol.GetBonds():
            start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bond_type_idx = BOND_INDEX[bond.GetBondType()]
            bond_dir_idx = BONDDIR_INDEX[bond.GetBondDir()]

            edge_index0[edge_offset] = start
            edge_index1[edge_offset] = end
            edge_attr0[edge_offset] = bond_type_idx
            edge_attr1[edge_offset] = bond_dir_idx
            edge_offset += 1

            edge_index0[edge_offset] = end
            edge_index1[edge_offset] = start
            edge_attr0[edge_offset] = bond_type_idx
            edge_attr1[edge_offset] = bond_dir_idx
            edge_offset += 1
        bond_feature_end = time.perf_counter() if self.enable_profiling else None

        if bond_count > 0:
            edge_index = torch.from_numpy(edge_index_np)
            edge_attr = torch.from_numpy(edge_attr_np)

        profile_timings = None
        if self.enable_profiling:
            profile_timings = (
                parse_end - parse_start,
                add_hs_end - add_hs_start,
                atom_feature_end - atom_feature_start,
                bond_feature_end - bond_feature_start,
            )

        return {
            "x": x,
            "edge_index": edge_index,
            "edge_attr": edge_attr,
            "has_unknown_atom": has_unknown_atom,
        }, profile_timings

    def __getitem__(self, index):
        total_start = time.perf_counter() if self.enable_profiling else None
        raw_smiles = self.smiles_data[index]
        graph_payload = self._graph_cache.get(raw_smiles)
        profile_timings = None

        if graph_payload is None:
            graph_payload, profile_timings = self._build_graph_payload(raw_smiles)
            self._graph_cache[raw_smiles] = graph_payload
        elif self.enable_profiling:
            profile_timings = (0.0, 0.0, 0.0, 0.0)

        if graph_payload["has_unknown_atom"]:
            print(self.id_data[index])

        data = Data(
            x=graph_payload["x"].clone(),
            edge_index=graph_payload["edge_index"].clone(),
            edge_attr=graph_payload["edge_attr"].clone(),
            chem_id=self.id_data[index],
        )
        if self.enable_profiling:
            data.graph_profile = torch.tensor(
                [[
                    profile_timings[0],
                    profile_timings[1],
                    profile_timings[2],
                    profile_timings[3],
                    time.perf_counter() - total_start,
                ]],
                dtype=torch.float32,
            )
        
        return data

    def __len__(self):
        return len(self.smiles_data)
    
    def get(self, index):
        return self.__getitem__(index)

    def len(self):
        return self.__len__()
    

# Function to generate the molecular scaffolds
def _generate_scaffold(smiles, include_chirality=False):
    mol = Chem.MolFromSmiles(smiles)
    scaffold = MurckoScaffoldSmiles(mol=mol, includeChirality=include_chirality)
    return scaffold


# Function to separate structures based on scaffolds
def generate_scaffolds(smile_list):

    """
    Generate molecular MURCKO scaffolds from a list of SMILES strings.

    Parameters:
    - smile_list (list): List of SMILES strings.

    Returns:
    - scaffold_sets (list): List of scaffold sets.
    """

    scaffolds = {}
    data_len = len(smile_list)

    print("About to generate scaffolds")
    for ind, smiles in enumerate(smile_list):
        scaffold = _generate_scaffold(smiles)
        if scaffold not in scaffolds:
            scaffolds[scaffold] = [ind]
        else:
            scaffolds[scaffold].append(ind)

    # Sort from largest to smallest scaffold sets
    scaffolds = {key: sorted(value) for key, value in scaffolds.items()}
    scaffold_sets = [
        scaffold_set for (scaffold, scaffold_set) in sorted(
            scaffolds.items(), key=lambda x: (len(x[1]), x[1][0]), reverse=True)
    ]
    return scaffold_sets

# Separate train, validation and test sets based on scaffolds
def scaffold_split(data_df, valid_size, test_size, smile_column, id_column):
    """
    Split data based on molecular scaffolds.

    Parameters:
    - data_df (pandas.DataFrame): DataFrame containing data to split.
    - valid_size (float): Proportion of data to allocate for validation.
    - test_size (float): Proportion of data to allocate for testing.
    - smile_column (str): Name of the column containing SMILES strings.
    - id_column (str): Name of the column containing molecule IDs.

    Returns:
    - train_ids (list): List of molecule IDs for the training set.
    - valid_ids (list): List of molecule IDs for the validation set.
    - test_ids (list): List of molecule IDs for the test set.
    """

    # Determine molecular scaffolds
    dataset = data_df[smile_column].tolist()
    train_size = 1.0 - valid_size - test_size
    scaffold_sets = generate_scaffolds(dataset)

    # Determine splits
    train_cutoff = train_size * len(dataset)
    valid_cutoff = (train_size + valid_size) * len(dataset)
    train_inds: List[int] = []
    valid_inds: List[int] = []
    test_inds: List[int] = []

    print("About to sort in scaffold sets")
    for scaffold_set in scaffold_sets:
        if len(train_inds) + len(scaffold_set) > train_cutoff:
            if len(train_inds) + len(valid_inds) + len(scaffold_set) > valid_cutoff:
                test_inds += scaffold_set
            else:
                valid_inds += scaffold_set
        else:
            train_inds += scaffold_set

    # Gather chem_ids based on 
    chemical_ids = data_df[id_column].tolist()
    train_ids = [chemical_ids[ind] for ind in train_inds]
    valid_ids = [chemical_ids[ind] for ind in valid_inds]
    test_ids = [chemical_ids[ind] for ind in test_inds]

    return train_ids, valid_ids, test_ids
    

# Function to split the dataset
def split_dataset(smile_df, valid_size, test_size, split_strategy, smile_col, id_col):

    """
    Split dataset into training, validation, and test sets.

    Parameters:
    - smile_df (pandas.DataFrame): DataFrame containing SMILES data.
    - valid_size (float): Proportion of data to allocate for validation.
    - test_size (float): Proportion of data to allocate for testing.
    - split_strategy (str): Splitting strategy ("random" or "scaffold").
    - smile_col (str): Name of the column containing SMILES strings.
    - id_col (str): Name of the column containing molecule IDs.

    Returns:
    - splitted_smiles_df (pandas.DataFrame): DataFrame with split information.
    """

    # Determine the splitting strategy
    if split_strategy == "random":

        # Determine the number of samples
        n_samples = smile_df.shape[0]

        # Randomly shuffle the indices
        chem_ids = smile_df[id_col].values
        np.random.shuffle(chem_ids)

        # Grab the validation ids
        valid_split = int(np.floor(valid_size * n_samples))
        valid_ids = chem_ids[:valid_split]

        # Grab the test ids
        test_split = int(np.floor(test_size * n_samples))
        test_ids = chem_ids[valid_split:(valid_split + test_split)]

        # Grab the train ids
        train_ids = chem_ids[(valid_split + test_split):]

    elif split_strategy == "scaffold":
        train_ids, valid_ids, test_ids = scaffold_split(smile_df, valid_size, test_size, smile_column=smile_col, id_column=id_col)

    # Add column with split information
    smile_df["split"]  = smile_df[id_col].apply(lambda x: "train" if x in train_ids else "valid" if x in valid_ids else "test")

    return smile_df

# Function to generate the molecular representation with MolE
def _resolve_graph_worker_count(num_graph_workers, dataset_size):
    if dataset_size <= 0:
        return 0

    if isinstance(num_graph_workers, str):
        if num_graph_workers != "auto":
            num_graph_workers = int(num_graph_workers)
        else:
            available_cpus = os.cpu_count() or 1
            num_graph_workers = max(1, available_cpus // 2)
            if dataset_size <= 1:
                num_graph_workers = 0

    if num_graph_workers is None:
        return 0

    return max(0, min(int(num_graph_workers), dataset_size))


def _resolve_prefetch_batches(prefetch_batches, num_graph_workers):
    if num_graph_workers <= 0:
        return None

    if isinstance(prefetch_batches, str):
        if prefetch_batches != "auto":
            prefetch_batches = int(prefetch_batches)
        else:
            prefetch_batches = 2

    if prefetch_batches is None:
        return 2

    return max(1, int(prefetch_batches))


def _is_cuda_device(device):
    if isinstance(device, torch.device):
        return device.type == "cuda"

    return str(device).startswith("cuda")


@contextmanager
def _representation_determinism(device, enabled=False):
    if not enabled or not _is_cuda_device(device):
        yield
        return

    with _DETERMINISTIC_TORCH_LOCK:
        original_algorithms = torch.are_deterministic_algorithms_enabled()
        original_cudnn_deterministic = torch.backends.cudnn.deterministic
        original_cudnn_benchmark = torch.backends.cudnn.benchmark

        torch.use_deterministic_algorithms(True)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

        try:
            yield
        finally:
            torch.use_deterministic_algorithms(original_algorithms)
            torch.backends.cudnn.deterministic = original_cudnn_deterministic
            torch.backends.cudnn.benchmark = original_cudnn_benchmark


_REPRESENTATION_PROFILE_SUM_FIELDS = (
    "graph_items",
    "rdkit_parse_seconds",
    "add_hs_seconds",
    "atom_feature_seconds",
    "bond_feature_seconds",
    "graph_total_seconds",
    "dataloader_setup_seconds",
    "dataloader_iter_seconds",
    "model_forward_seconds",
)


def _new_representation_profile(graph_batch_size, graph_workers):
    return {
        "graph_items": 0,
        "rdkit_parse_seconds": 0.0,
        "add_hs_seconds": 0.0,
        "atom_feature_seconds": 0.0,
        "bond_feature_seconds": 0.0,
        "graph_total_seconds": 0.0,
        "dataloader_setup_seconds": 0.0,
        "dataloader_iter_seconds": 0.0,
        "model_forward_seconds": 0.0,
        "graph_batch_size": graph_batch_size,
        "graph_workers": graph_workers,
    }


def _accumulate_representation_profile(total_profile, batch_profile):
    for field in _REPRESENTATION_PROFILE_SUM_FIELDS:
        total_profile[field] += batch_profile.get(field, 0.0)

    if "graph_batch_size" in batch_profile:
        total_profile["graph_batch_size"] = batch_profile["graph_batch_size"]
    if "graph_workers" in batch_profile:
        total_profile["graph_workers"] = batch_profile["graph_workers"]


def _build_representation_loader(
    smile_df,
    column_str,
    id_str,
    batch_size,
    num_graph_workers,
    prefetch_batches,
    device,
    enable_profiling,
):
    dataloader_setup_start = time.perf_counter() if enable_profiling else None
    molecular_graph_dataset = MoleculeDataset(
        smile_df,
        column_str,
        id_str,
        enable_profiling=enable_profiling,
    )
    dataset_size = len(molecular_graph_dataset)
    graph_batch_size = int(batch_size)
    resolved_workers = _resolve_graph_worker_count(num_graph_workers, dataset_size)
    resolved_prefetch = _resolve_prefetch_batches(prefetch_batches, resolved_workers)
    use_cuda_transfer = _is_cuda_device(device)

    loader_kwargs = {
        "dataset": molecular_graph_dataset,
        "batch_size": graph_batch_size,
        "shuffle": False,
        "num_workers": resolved_workers,
    }
    if use_cuda_transfer:
        loader_kwargs["pin_memory"] = True
    if resolved_prefetch is not None:
        loader_kwargs["prefetch_factor"] = resolved_prefetch
    if resolved_workers > 0:
        loader_kwargs["persistent_workers"] = True

    dataloader_setup_seconds = 0.0
    if enable_profiling:
        dataloader_setup_seconds = time.perf_counter() - dataloader_setup_start

    return (
        DataLoader(**loader_kwargs),
        graph_batch_size,
        resolved_workers,
        use_cuda_transfer,
        dataloader_setup_seconds,
    )


def _profile_graph_batch(batch_profile, graph_batch):
    graph_profile = getattr(graph_batch, "graph_profile", None)
    if graph_profile is None:
        return

    graph_profile = graph_profile.view(-1, 5).sum(dim=0).tolist()
    batch_profile["graph_items"] += int(len(getattr(graph_batch, "chem_id", [])))
    batch_profile["rdkit_parse_seconds"] += float(graph_profile[0])
    batch_profile["add_hs_seconds"] += float(graph_profile[1])
    batch_profile["atom_feature_seconds"] += float(graph_profile[2])
    batch_profile["bond_feature_seconds"] += float(graph_profile[3])
    batch_profile["graph_total_seconds"] += float(graph_profile[4])


def iter_batch_representation(
    smile_df,
    dl_model,
    column_str,
    id_str,
    batch_size=10_000,
    id_is_str=True,
    device="cuda:0",
    num_graph_workers=0,
    graph_batch_size=None,
    prefetch_batches=2,
    enable_profiling=False,
    deterministic_representation=False,
):

    """
    Yield molecular representation mini-batches with optional profiling.

    Each yielded item contains:
    - batch_index: incremental mini-batch index
    - chem_ids: list of chem_ids in the yielded embedding batch
    - embedding_batch: pandas.DataFrame with embeddings indexed by chem_id
    - profiling: per-batch profiling summary when enabled, otherwise None
    """

    del id_is_str

    resolved_batch_size = int(graph_batch_size or batch_size)
    (
        molecular_graph_loader,
        resolved_batch_size,
        resolved_workers,
        use_cuda_transfer,
        dataloader_setup_seconds,
    ) = _build_representation_loader(
        smile_df=smile_df,
        column_str=column_str,
        id_str=id_str,
        batch_size=resolved_batch_size,
        num_graph_workers=num_graph_workers,
        prefetch_batches=prefetch_batches,
        device=device,
        enable_profiling=enable_profiling,
    )

    dl_model.eval()
    with _representation_determinism(device, deterministic_representation), torch.no_grad():
        loader_iter = iter(molecular_graph_loader)
        batch_index = 0
        while True:
            batch_profile = None
            iter_start = time.perf_counter() if enable_profiling else None
            try:
                graph_batch = next(loader_iter)
            except StopIteration:
                break

            if enable_profiling:
                batch_profile = _new_representation_profile(
                    graph_batch_size=resolved_batch_size,
                    graph_workers=resolved_workers,
                )
                batch_profile["dataloader_setup_seconds"] = (
                    dataloader_setup_seconds if batch_index == 0 else 0.0
                )
                batch_profile["dataloader_iter_seconds"] = time.perf_counter() - iter_start
                _profile_graph_batch(batch_profile, graph_batch)

            forward_start = time.perf_counter() if enable_profiling else None
            graph_batch = graph_batch.to(device, non_blocking=use_cuda_transfer)
            h_representation, _ = dl_model(graph_batch)
            if enable_profiling:
                batch_profile["model_forward_seconds"] = time.perf_counter() - forward_start

            chem_ids = list(graph_batch.chem_id)
            batch_df = pd.DataFrame(h_representation.cpu().numpy(), index=chem_ids)
            yield {
                "batch_index": batch_index,
                "chem_ids": chem_ids,
                "embedding_batch": batch_df,
                "profiling": batch_profile,
            }
            batch_index += 1


def batch_representation(
    smile_df,
    dl_model,
    column_str,
    id_str,
    batch_size=10_000,
    id_is_str=True,
    device="cuda:0",
    num_graph_workers=0,
    graph_batch_size=None,
    prefetch_batches=2,
    enable_profiling=False,
    deterministic_representation=False,
):

    """
    Generate molecular representations using a Deep Learning model.

    Parameters:
    - smile_df (pandas.DataFrame): DataFrame containing SMILES data.
    - dl_model: Deep Learning model for molecular representation.
    - column_str (str): Name of the column containing SMILES strings.
    - id_str (str): Name of the column containing molecule IDs.
    - batch_size (int, optional): Batch size for processing (default is 10,000).
    - id_is_str (bool, optional): Whether IDs are strings (default is True).
    - device (str, optional): Device for computation (default is "cuda:0").

    Returns:
    - chem_representation (pandas.DataFrame): DataFrame containing molecular representations.
    """
    
    resolved_batch_size = int(graph_batch_size or batch_size)
    profiling = (
        _new_representation_profile(
            graph_batch_size=resolved_batch_size,
            graph_workers=0,
        )
        if enable_profiling
        else None
    )
    batch_dataframes = []
    for batch in iter_batch_representation(
        smile_df=smile_df,
        dl_model=dl_model,
        column_str=column_str,
        id_str=id_str,
        batch_size=batch_size,
        id_is_str=id_is_str,
        device=device,
        num_graph_workers=num_graph_workers,
        graph_batch_size=graph_batch_size,
        prefetch_batches=prefetch_batches,
        enable_profiling=enable_profiling,
        deterministic_representation=deterministic_representation,
    ):
        batch_dataframes.append(batch["embedding_batch"])
        if enable_profiling and profiling is not None and batch["profiling"] is not None:
            _accumulate_representation_profile(profiling, batch["profiling"])

    # Concatenate the dataframes
    chem_representation = pd.concat(batch_dataframes)
    if enable_profiling and profiling is not None:
        chem_representation.attrs["profiling"] = profiling

    return chem_representation

# Function to load a pre-trained model
def load_pretrained_model(pretrain_architecture, pretrained_model, pretrained_dir = "../pretrained_model", device="cuda:0"):

    """
    Load a pre-trained MolE model.

    Parameters:
    - pretrain_architecture (str): Architecture of the pre-trained model.
    - pretrained_model (str): Name of the pre-trained MolE model.
    - pretrained_dir (str, optional): Directory containing pre-trained models (default is "../pretrained_model").
    - device (str, optional): Device for computation (default is "cuda:0").

    Returns:
    - model: Loaded pre-trained model.
    """

    # Read model configuration
    config = yaml.load(open(os.path.join(pretrained_dir, pretrained_model, "config.yaml"), "r"), Loader=yaml.FullLoader)
    model_config = config["model"]

    # Instantiate model
    if pretrain_architecture == "gin_concat":
        from models.ginet_concat import GINet
        model = GINet(**model_config).to(device)
    
    # Load pre-trained weights
    model_pth_path = os.path.join(pretrained_dir, pretrained_model, "model.pth")
    print(model_pth_path)

    state_dict = torch.load(model_pth_path, map_location=device)
    model.load_my_state_dict(state_dict)

    return model

# Function to generate the ECFP4 as an array
def fp_array(fingerprin_object):

    """
    Convert fingerprint object to NumPy array.

    Parameters:
    - fingerprin_object: Fingerprint object.

    Returns:
    - array (numpy.ndarray): NumPy array representation of the fingerprint.
    """

    # Initialise an array full of zeros
    array = np.zeros((0,), dtype=np.int8)

    # Dump fingerprint info into array
    DataStructs.ConvertToNumpyArray(fingerprin_object, array)

    return array

# Function to generate the ECFP4 representation
def generate_fps(smile_df, smile_col, id_col):

    """
    Generate Extended-Connectivity Fingerprints (ECFP4) representations for molecules.

    Parameters:
    - smile_df (pandas.DataFrame): DataFrame containing SMILES data.
    - smile_col (str): Name of the column containing SMILES strings.
    - id_col (str): Name of the column containing molecule IDs.

    Returns:
    - fps_dataframe (pandas.DataFrame): DataFrame containing ECFP4 representations.
    """

    # Generate fingerprints
    mol_objs = [Chem.MolFromSmiles(smile) for smile in smile_df[smile_col].tolist()]
    fp_objs = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in mol_objs]

    # Place fingerprints in array
    fps_arrays = [fp_array(fp) for fp in fp_objs]

    # Create dataframe
    fps_matrix = np.stack(fps_arrays, axis=0 )
    fps_dataframe = pd.DataFrame(fps_matrix, index=smile_df[id_col].tolist())

    return fps_dataframe

def generate_descriptors(smile_df, smile_col, id_col):

    """
    Generate molecular descriptors for molecules.

    Parameters:
    - smile_df (pandas.DataFrame): DataFrame containing SMILES data.
    - smile_col (str): Name of the column containing SMILES strings.
    - id_col (str): Name of the column containing molecule IDs.

    Returns:
    - chemdesc_df (pandas.DataFrame): DataFrame containing molecular descriptors.
    """


    # Iterate over chemicals, making sure to catch exceptions
    descriptor_list = []
    for smile, id in zip(smile_df[smile_col].tolist(), smile_df[id_col].tolist()):
        try:
             desc_df = calc_descriptors(smile, id)
             descriptor_list.append(desc_df)

        except:
            print(f"Could not compute descriptors for {id}")
        
        
    # Concatenate and output
    chemdesc_df = pd.concat([d for d in descriptor_list if type(d) != str])
    chemdesc_df = chemdesc_df.rename(columns = {"chem_id": id_col}).set_index(id_col)

    return chemdesc_df

# Main function to process the dataset
def process_dataset(dataset_path,
                    smile_column_str = "rdkit_no_salt", 
                    id_column_str = "prestwick_ID",
                    pretrain_architecture=None, 
                    pretrained_model = None, 
                    split_approach="scaffold", 
                    validation_proportion=0.1, 
                    test_proportion=0.1, 
                    dataset_split=True,
                    device="cuda:0"):
    
    """
    Process the dataset to generate molecular representations.

    Parameters:
    - dataset_path (str): Path to the dataset file.
    - pretrain_architecture (str): Architecture of the pre-trained model or method ("ECFP4", "ChemDesc", or custom).
    - pretrained_model (str): Name of the pre-trained model. Can also be "MolCLR" or "ECFP4".
    - split_approach (str, optional): Splitting approach ("scaffold" or "random") (default is "scaffold").
    - validation_proportion (float, optional): Proportion of data to allocate for validation (default is 0.1).
    - test_proportion (float, optional): Proportion of data to allocate for testing (default is 0.1).
    - smile_column_str (str, optional): Name of the column containing SMILES strings (default is "rdkit_no_salt").
    - id_column_str (str, optional): Name of the column containing molecule IDs (default is "prestwick_ID").
    - dataset_split (bool, optional): Whether to split the dataset into train, validation, and test sets (default is True).
    - device (str): Device to use for computation (default is "cuda:0"). Can also be "cpu".

    Returns:
    - splitted_smiles_df (pandas.DataFrame): DataFrame with split information if split_data=True.
    - udl_representation (pandas.DataFrame): DataFrame containing molecular representations if split_data=False.
    """

    # First we read in the smiles as a dataframe
    smiles_df = read_smiles(dataset_path, smile_col=smile_column_str, id_col=id_column_str)

    # The we split the dataset into train, validation and test
    if dataset_split:
        splitted_smiles_df = split_dataset(smiles_df, validation_proportion, test_proportion, split_approach, smile_column_str, id_column_str)

    # Determine the representation
    if pretrained_model == "ECFP4":
        udl_representation = generate_fps(smiles_df, smile_column_str, id_column_str)
        
    elif pretrained_model == "ChemDesc":
        udl_representation = generate_descriptors(smiles_df, smile_column_str, id_column_str)
    
    else:
        # Now we load our pretrained model
        pmodel = load_pretrained_model(pretrain_architecture, pretrained_model, device=device)
        # Obtain the requested representation
        udl_representation = batch_representation(smiles_df, pmodel, smile_column_str, id_column_str, device=device)

    if dataset_split:
        return splitted_smiles_df, udl_representation
    
    else:
        return udl_representation
