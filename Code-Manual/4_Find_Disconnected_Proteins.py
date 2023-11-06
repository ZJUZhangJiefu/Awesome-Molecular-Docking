import glob
import os

import networkx as nx
from biopandas.pdb import PandasPdb

from scipy import spatial
from tqdm import tqdm
import numpy as np
import pandas as pd
def write_strings_to_txt(strings: list, path):
    # every string of the list will be saved in one line
    textfile = open(path, "w")
    for element in strings:
        textfile.write(element + "\n")
    textfile.close()
pdb_path = 'data/PDBBind'
pdbbind_names = os.listdir(pdb_path)
# read the the index file of PDBBind v2020
df_pdb_id = pd.read_csv('data/PDBbind_index/INDEX_general_PL_name.2020', sep="  ", comment='#', header=None, names=['complex_name', 'year', 'pdb_id', 'd', 'e','f','g','h','i','j','k','l','m','n','o'])
# select the attributes named ['complex_name','year','pdb_id']
# year is used in EquiBind, TankBind and DiffDock to split the train and test set of PDBBind v2020
df_pdb_id = df_pdb_id[['complex_name','year','pdb_id']]

# read and select some other data
df_data = pd.read_csv('data/PDBbind_index/INDEX_general_PL_data.2020', sep="  ", comment='#', header=None, names=['complex_name','resolution','year', 'logkd', 'kd', 'reference', 'ligand_name', 'a', 'b', 'c'])
df_data = df_data[['complex_name','resolution','year', 'logkd', 'kd', 'reference', 'ligand_name']]

# distance cutoff: 5A, used to judge whether a protain is connected
cutoff = 5
connected = []

for name in tqdm(pdbbind_names):
    # read the protein pdb files that have been processed by openbabel and reduce in file 1 and 2
    # and select their atoms  
    df = PandasPdb().read_pdb(os.path.join(pdb_path, name, f'{name}_protein_obabel_reduce.pdb')).df['ATOM']
    # rename the attributes for better understanding
    df.rename(columns={'chain_id': 'chain', 'residue_number': 'residue', 'residue_name': 'resname',
                       'x_coord': 'x', 'y_coord': 'y', 'z_coord': 'z', 'element_symbol': 'element'}, inplace=True)
    # group the protein atoms by different chains
    df = list(df.groupby(['chain']))  ## Not the same as sequence order !

    chain_coords_list = []
    for chain in df:
        # get protein chains
        chain_coords_list.append(chain[1][['x', 'y', 'z']].to_numpy().squeeze().astype(np.float32))

    num_chains = len(chain_coords_list)
    distance = np.full((num_chains, num_chains), -np.inf)
    # calculate the minimum distance between each pair of chains in the protein
    for i in range(num_chains - 1):
        for j in range((i + 1), num_chains):
            pairwise_dis = spatial.distance.cdist(chain_coords_list[i],chain_coords_list[j])
            distance[i, j] = np.min(pairwise_dis)
            distance[j, i] = np.min(pairwise_dis)
    src_list = []
    dst_list = []
    
    # select protein chains whose distance < cutoff
    # source chain -> destination chain
    for i in range(num_chains):
        dst = list(np.where(distance[i, :] < cutoff)[0])
        src = [i] * len(dst)
        src_list.extend(src)
        dst_list.extend(dst)
    # use information source chain -> destination chain to construct a graph
    # distance between each pair of source chain and destination chain < cutoff(5A)
    graph = nx.Graph()
    graph.add_edges_from(zip(src_list, dst_list))
    # the graph is connected means: 
    # do any division of the protein into 2 parts
    # (divide it by chains, some chains form a part, and the other chains form another part)
    # then the min distance of 2 parts is always smaller than cutoff(5A)
    if nx.is_connected(graph):
        connected.append(name)
    else:
        print(f'not connected: {name}')
# write the names of complexes (in which the protein is connected) into a new txt file
write_strings_to_txt(connected, f'data/complex_names_connected_by_{cutoff}')
print(len(connected))
