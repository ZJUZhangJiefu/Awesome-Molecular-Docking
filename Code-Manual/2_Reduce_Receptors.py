# See readme for details about installing "reduce"
import os
import subprocess
import time
from tqdm import tqdm
start_time = time.time()
data_path = 'data/PDBBind'
overwrite = False
names = sorted(os.listdir(data_path))

for i, name in tqdm(enumerate(names)):
    rec_path = os.path.join(data_path, name, f'{name}_protein_obabel.pdb')
    # reduce -Trim: remove (rather than add) hydrogens and skip all optimizations
    return_code = subprocess.run(
        f"reduce -Trim {rec_path} > {os.path.join(data_path, name, f'{name}_protein_obabel_reduce_tmp.pdb')}", shell=True)
    print(return_code)
    # reduce -HIS: create NH hydrogens on HIS rings
    return_code2 = subprocess.run(
        f"reduce -HIS {os.path.join(data_path, name, f'{name}_protein_obabel_reduce_tmp.pdb')} > {os.path.join(data_path, name, f'{name}_protein_obabel_reduce.pdb')}", shell=True)
    print(return_code2)
    # remove the temporary files produced when running the above code
    return_code2 = subprocess.run(
        f"rm {os.path.join(data_path, name, f'{name}_protein_obabel_reduce_tmp.pdb')}",
        shell=True)
    print(return_code2)

print("--- %s seconds ---" % (time.time() - start_time))
