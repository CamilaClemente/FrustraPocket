# FrustraPocket

The main goal of FrustraPocket is to predict protein-ligand and catalytic pockets and perform molecular docking of a specific ligand to each predicted pocket.
We highlight that the important aspect of this pipeline is that it detects pockets of a protein, using frustration, and calculates its center of mass to performed a guided docking. So, you don't need to know the coordinates of the pockets to be analyzed!

Tools required to predict protein-ligand and catalytic sites (1):
FrustratometerR: https://github.com/proteinphysiologylab/frustratometeR 

For visualization:
Pymol: https://pymol.org/ (sudo apt-get install pymol)
	
	
The input files are: a PDBId 

The pipeline is code in python3. 
You can run the pipeline as a follow:

$ python3 main_fp.py --JobId 1a0i --PdbId 1a0i

For a particular chain

$ python3 main_fp.py --JobId 1a0i --PdbId 1a0i --chain A

If you don't have the PDB structures the pipeline will download the PDB.
