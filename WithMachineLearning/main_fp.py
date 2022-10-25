import os
import sys
from os.path import exists
sys.path.append('/home/maria/Documentos/')#Path to Functions.py file
import Functions
import argparse

parser = argparse.ArgumentParser(description='Calculation of the frustration logo.')
parser.add_argument("--JobId", help="The name of the job")
parser.add_argument("--PdbId", help="The ID of the pdb")
args = parser.parse_args()

if not exists(args.PdbId+'.pdb'):
	os.system('wget https://files.rcsb.org/view/'+args.PdbId+'.pdb -O '+args.PdbId+'.pdb')

chains=Functions.fst_files(args.JobId,args.PdbId)

for x in chains:
	Functions.runmodel(args.JobId,args.PdbId,x)
	Functions.output(args.JobId,args.PdbId,x)

