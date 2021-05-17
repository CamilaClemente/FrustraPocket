# FrustraPocket

The main goal of FrustraPocket is to predict protein-ligand and catalytic pockets and perform molecular docking of a specific ligand to each predicted pocket.
We highlight that the important aspect of this pipeline is that it detects pockets of a protein, using frustration, and calculates its center of mass to performed a guided docking. So, you don't need to know the coordinates of the pockets to be analyzed!

You can use the pipeline in two ways:
1) Only for predict protein-ligand and catalityc pockets
2) To predict protein-ligand and catalytic pockets and run a docking with a specific ligand to each pocket predicted.

Tools required to predict protein-ligand and catalytic sites (1):
FrustratometerR: https://github.com/proteinphysiologylab/frustratometeR 

Additional tool to predict protein-ligand and catalytic sites and run a docking with a specific ligand to each pocket predicted (2):
MGL Tools: http://mgltools.scripps.edu/downloads) - 
Autodock Vina: http://vina.scripps.edu/download.html

For visualization:
Pymol: https://pymol.org/ (sudo apt-get install pymol)
	
	
The input files are: a PDBId and a ligand (PDB Format) to do the docking (only if it is necessary).
The output files:
	- The frustration calculations (mutational)
	- The pockets predicted and a pymol script to vizualize the pockets in the protein structure
	- The center of mass for each pocket predicted
	- If you run a docking, the vina results for each pocket predicted


The pipeline is code in python3. 
You can run the pipeline as a follow:
Only for pocket prediction:

$ python3 RunFrustraPocket.py PDBId (python3 RunFrustraPocket.py 4new)

For pocket prediction and docking (the ligand in PDB format must be in the folder where the RunFrustraPocket.py is located):

$ python3 RunFrustraPocket.py PDBId Ligand (python3 RunFrustraPocket.py 4new FAD) 

You have to change a path in the RunFrustraPocket.py script:
- line 203, pythonpath="pythonsh_path_from_MGLTools"


If you have a dataset of the proteins to run in PDB format and you want to use them, you have to define a PATH where the structures are, to do this in RunFrustraPocket.py in line 119 set the pathPDB="path_to_PDB_Files"

If you don't have the PDB structures the pipeline will download the PDB.
