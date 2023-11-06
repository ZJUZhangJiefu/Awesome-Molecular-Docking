"""
First, use openbabel to process the receptors.
"""
import os
import subprocess
import time
from tqdm import tqdm
# record the start time 
start_time = time.time()
# set your own data path of PDBBind
data_path = 'data/PDBBind'
"""
    Your directory of PDBBind dataset should be like this:  
    -PDBBind 
      -1a0t
        -1a0t_ligand.mol2
        -1a0t_ligand.sdf
        -1a0t_pocket.pdb
        -1a0t_protein.pdb
      -1a1c
      ...
"""
overwrite = False
# obtain the name of protein-ligand complexes, such as 1a0t, etc.
names = sorted(os.listdir(data_path))
for i, name in tqdm(enumerate(names)):
    # the path of receptor(small drug molecule)
    rec_path = os.path.join(data_path, name, f'{name}_protein.pdb')
    # process the receptor file whose path is rec_path, and output it to f'{name}_protein_obabel.pdb
    return_code = subprocess.run(
        f"obabel {rec_path} -O{os.path.join(data_path, name, f'{name}_protein_obabel.pdb')}", shell=True)
    print(return_code)
