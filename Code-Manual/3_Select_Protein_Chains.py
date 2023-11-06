import os
import warnings
import numpy as np
import prody
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from scipy import spatial
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
def read_molecule(molecule_file, sanitize=False, calc_charges=False, remove_hs=False):
    """Load a molecule from a file of format ``.mol2`` or ``.sdf`` or ``.pdbqt`` or ``.pdb``.

    Parameters
    ----------
    molecule_file : str
        Path to file for storing a molecule, which can be of format ``.mol2`` or ``.sdf``
        or ``.pdbqt`` or ``.pdb``.
    sanitize : bool
        Whether sanitization is performed in initializing RDKit molecule instances. See
        https://www.rdkit.org/docs/RDKit_Book.html for details of the sanitization.
        Default to False.
    calc_charges : bool
        Whether to add Gasteiger charges via RDKit. Setting this to be True will enforce
        ``sanitize`` to be True. Default to False.
    remove_hs : bool
        Whether to remove hydrogens via RDKit. Note that removing hydrogens can be quite
        slow for large molecules. Default to False.
    use_conformation : bool
        Whether we need to extract molecular conformation from proteins and ligands.
        Default to True.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        RDKit molecule instance for the loaded molecule.
    coordinates : np.ndarray of shape (N, 3) or None
        The 3D coordinates of atoms in the molecule. N for the number of atoms in
        the molecule. None will be returned if ``use_conformation`` is False or
        we failed to get conformation information.
    """
    # create a mol object from different formats of molecule files
    if molecule_file.endswith('.mol2'):
        mol = Chem.MolFromMol2File(molecule_file, sanitize=False, removeHs=False)
    elif molecule_file.endswith('.sdf'):
        supplier = Chem.SDMolSupplier(molecule_file, sanitize=False, removeHs=False)
        mol = supplier[0]
    elif molecule_file.endswith('.pdbqt'):
        with open(molecule_file) as file:
            pdbqt_data = file.readlines()
        pdb_block = ''
        for line in pdbqt_data:
            pdb_block += '{}\n'.format(line[:66])
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    elif molecule_file.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(molecule_file, sanitize=False, removeHs=False)
    else:
        return ValueError('Expect the format of the molecule_file to be '
                          'one of .mol2, .sdf, .pdbqt and .pdb, got {}'.format(molecule_file))

    try:
        if sanitize or calc_charges:
            # sanitize the molecule
            Chem.SanitizeMol(mol)

        if calc_charges:
            # Compute Gasteiger charges on the molecule.
            try:
                AllChem.ComputeGasteigerCharges(mol)
            except:
                warnings.warn('Unable to compute charges for the molecule.')

        if remove_hs:
            # remove Hydrogen atoms on the molecule
            mol = Chem.RemoveHs(mol, sanitize=sanitize)
    except:
        return None

    return mol

cutoff = 10
data_dir = 'data/PDBBind'
names = os.listdir(data_dir)

io = PDBIO()
biopython_parser = PDBParser()

for name in tqdm(names):
    # recepter path
    # the postfix "_obabel_reduce" means:
    # the protein receptors have been preprocessed by openbabel and reduce
    rec_path = os.path.join(data_dir, name, f'{name}_protein_obabel_reduce.pdb')
    # read ligand molecule into a mol object from different file formats 
    lig = read_molecule(os.path.join(data_dir, name, f'{name}_ligand.sdf'), sanitize=True, remove_hs=False)
    if lig == None:
        lig = read_molecule(os.path.join(data_dir, name, f'{name}_ligand.mol2'), sanitize=True, remove_hs=False)
    if lig == None:
        print('ligand was none for ', name)
        with open('select_chains.log', 'a') as file:
            file.write(f'{name}\n')
        continue
    # the conformer of ligand
    conf = lig.GetConformer()
    # the 3D coordinates of each atom in the ligand
    # Note:
    # The ligand is represented atom-wise(fine-grained)
    # The receptor is represented residue-wise(coarse-grained)
    lig_coords = conf.GetPositions()
    # filter PDBConstructionWarning
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=PDBConstructionWarning)
        # parse the receptor
        structure = biopython_parser.get_structure('random_id', rec_path)
        rec = structure[0]
    
    min_distances = []
    coords = []
    valid_chain_ids = []
    lengths = []
    for i, chain in enumerate(rec):
        chain_coords = []  # num_residues, num_atoms, 3
        chain_is_water = False
        count = 0
        invalid_res_ids = []
        for res_idx, residue in enumerate(chain):
            if residue.get_resname() == 'HOH':
                chain_is_water = True
            residue_coords = []
            c_alpha, n, c = False, False, False
            # judge the type of atoms
            for atom in residue:
                if atom.name == 'CA':
                    c_alpha = True
                if atom.name == 'N':
                    n = True
                if atom.name == 'C':
                    c = True
                # Gather the coordinates of atoms in a protein residue
                residue_coords.append(list(atom.get_vector()))
            # only append residue if it is an amino acid 
            # and not some weired molecule that is part of the complex
            # The condition of being an amino acid:
            # Containing C_alpha, N and C Atoms at the same time
            if c_alpha and n and c:  
                chain_coords.append(np.array(residue_coords))
                count += 1
            # else:record the invalid protein chains 
            else:
                invalid_res_ids.append(residue.get_id())
        # detach the invalid residues
        # (some invalid residue may be H2O or other things)
        for res_id in invalid_res_ids:
            chain.detach_child(res_id)
        if len(chain_coords) > 0:
            all_chain_coords = np.concatenate(chain_coords, axis=0)
            # calculate pairwise distance between ligand atoms and protain atoms
            # see the above code: atom coordinate->residue_coords->chain_coords->all_chain_coords
            distances = spatial.distance.cdist(lig_coords, all_chain_coords)
            min_distance = distances.min()
        else:
            min_distance = np.inf
        if chain_is_water:
            min_distances.append(np.inf)
        else:
            min_distances.append(min_distance)
        lengths.append(count)
        coords.append(chain_coords)
        # min distance is within cutoff(default: 10A) and the chain is not a water chain
        if min_distance < cutoff and not chain_is_water:
            valid_chain_ids.append(chain.get_id())
    min_distances = np.array(min_distances)
    if len(valid_chain_ids) == 0:
        # no valid chains
        valid_chain_ids.append(np.argmin(min_distances))
    valid_coords = []
    valid_lengths = []
    invalid_chain_ids = []
    # the for loop may be unneeded here
    for i, chain in enumerate(rec):
        if chain.get_id() in valid_chain_ids:
            valid_coords.append(coords[i])
            valid_lengths.append(lengths[i])
        else:
            invalid_chain_ids.append(chain.get_id())
    # parse the current protein receptor(PDB file format)
    prot = prody.parsePDB(rec_path)
    # select the valid chains of proteins according to valid_chain_ids
    # each protein chain is identified by a chain_id
    try:
        sel = prot.select(' or '.join(map(lambda c: f'chain {c}', valid_chain_ids)))
    except:
        print("Problematic valid_chain_ids: ", valid_chain_ids)
        print("Encountered illegal chains, and keep all chains in this protein.")
        sel = prot
    """
    NOTICE: 
    When encountering chains that have no name, I skip selecting the corresponding protein chains, 
    and keep the whole protein for simplicity.
    You may need to modify the try-except code section for more professional reasons.
    """
    # write the protein with valid chain selected into another PDB file
    prody.writePDB(os.path.join(data_dir,name,f'{name}_protein_processed2.pdb'),sel)
