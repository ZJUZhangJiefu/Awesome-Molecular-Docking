<table><tr><td bgcolor= #E0FFF#>
<font face = "Times New Roman">
<font size=3>
  
#### Introduction
This is a repository of **molecular docking from beginner to expert**. It starts from **manuals of environment setup and basic operations to clean the PDBBind v2020 dataset**, and will contain many **tricky parts involved in molecular docking**, for example, **IEGMN**.  
Both **ipynb files with detailed comments and segmented sample output** and **python files to run directly** will be provided.  

Hope I will help you gain an insight to molecular docking to the best of my ability, and many thanks for your welcome.  

#### Setting up the environment 
Recommended environment: Visual Studio Code with Linux Ubuntu 22.04. 
##### 1 Basic Environment Setup  
First, run the following commands to set up an appropriate environment.  
Notice: Some high versions of torch-geometric, torch-scatter or relevant libraries may encounter an error like "has no object torch_csc_tensor". Line 4, 5 will install some lower versions of geometric libraries, and they will not bring similar problems.   
```
conda create --name diffdock python=3.9
conda activate diffdock
conda install pytorch==1.11.0 pytorch-cuda=11.7 -c pytorch -c nvidia
pip install torch-scatter==2.0.9 torch-sparse==0.6.13 torch-cluster==1.6.0 torch-spline-conv==1.2.1 torch-geometric==2.0.4
python -m pip install PyYAML scipy "networkx[default]" biopython rdkit-pypi e3nn spyrmsd pandas biopandas prody dgl wandb
```
Then, run the following code for protein sequence embedding and protein structure generation.  
```
pip install "fair-esm[esmfold]"
pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'
```
##### 2 Install Openbabel  
Run the following commands to install openbabel:  
```
conda install openbabel -c conda-forge
```  
Then, run ```obabel -V``` to confirm that you've successfully installed `obabel`.  
##### 3 The esm pytorch parameter file  
Download "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt" to /home/ubuntu/.cache/torch/hub/checkpoints/esm2_t33_650M_UR50D.pt.  
Here, `ubuntu` should be replaced by your own user name.  
##### 4 Install reduce  
First clone the codes:  
```
sudo git clone https://github.com/rlabduke/reduce.git
```
Then, open the directory `reduce` and run `make` command.   
Finally, run `sudo make install` command in directory `reduce`.  
