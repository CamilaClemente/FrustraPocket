#!/usr/bin/python3
import os
import sys
from numpy import loadtxt
from tensorflow.keras.models import load_model
import pandas as pd
import pickle
from sklearn.metrics import accuracy_score, precision_score, recall_score, cohen_kappa_score, f1_score, confusion_matrix, roc_auc_score, roc_curve, auc, mean_absolute_error, mean_squared_error
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import xgboost as xgb

def mkdir(jobID):#make the necessary folders
	if os.path.isdir(jobID):
		os.system('rm -r '+jobID)
	os.system('mkdir '+jobID)
	os.system('mkdir '+jobID+'/Frustration')
	os.system('mkdir '+jobID+'/Frustration/FrustrationResults')
	os.system('mkdir '+jobID+'/Data')
	os.system('mkdir '+jobID+'/OutPutFiles')
	
def split_pdb(pdbId,jobID,chains=[]):#split the pdbs into chains and return a vector with the chains
	pdb=open(pdbId+'.pdb')
	if len(chains) == 0:
		for line in pdb.readlines():
			if line[0:4] == 'ATOM':
				out=open(jobID+'/Frustration/'+pdbId+'_'+line[21]+'.pdb','a')
				out.write(line)
				out.close()
				if line[21] not in chains:
					chains.append(line[21])
	else:
		for line in pdb.readlines():
			if line[0:4] == 'ATOM' and line[21] == chains:
				out=open(jobID+'/Frustration/'+pdbId+'_'+line[21]+'.pdb','a')
				out.write(line)
				out.close()
	pdb.close()
	return chains
	
def Frustra_PDBs(jobID):# frustration calculation of thr pdbs
	frustra=open(jobID+'/Frustra.r','w')
	directory=os.getcwd()+'/'+jobID+'/'
	frustra.write('library(frustratometeR)\nPdbsDir <- \"'+directory+'Frustration/\"\nResultsDir <- \"'+directory+'Frustration/FrustrationResults/\"\ndir_frustration(PdbsDir = PdbsDir, Mode = \"mutational\", ResultsDir = ResultsDir)')
	frustra.close()
	os.system('cd '+jobID+';Rscript Frustra.r')

def fst_files(jobID,pdb,chain=[]):
	mkdir(jobID)
	if len(chain) ==0:
		chains=split_pdb(pdb,jobID)
	else:
		chains=split_pdb(pdb,jobID,chain)
	Frustra_PDBs(jobID)
	for x in chains:
		#jobID+'/Frustration/'+pdb+'_'+x+'.pdb' /home/maria/Documentos/FrustraPocket/1a0i/Frustration/FrustrationResults/1a0i_A.done/FrustrationData/1a0i_A.pdb_mutational_5adens
		fst=open(jobID+'/Frustration/FrustrationResults/'+pdb+'_'+x+'.done/FrustrationData/'+pdb+'_'+x+'.pdb_mutational_5adens','r')
		resn=[]
		nocont=[]
		pohf=[]
		ponf=[]
		pomf=[]
		for line in fst.readlines():
			line=line.rstrip('\n')
			spfst=line.split()
			if spfst[0] != 'Res':
				resn.append(spfst[0])
				nocont.append(spfst[2])
				pohf.append(spfst[6])
				ponf.append(spfst[7])
				pomf.append(spfst[8])
		fst.close()
		os.system('cd '+jobID+'/Frustration/FrustrationResults/'+pdb+'_'+x+'.done/FrustrationData/;awk \'{print $1,$5}\' '+pdb+'_'+x+'.pdb_mutational > LD_aux; awk \'{print $2,$6}\' '+pdb+'_'+x+'.pdb_mutational >> LD_aux; sort -u -k1n LD_aux > LD;rm LD_aux')
		locald=open(jobID+'/Frustration/FrustrationResults/'+pdb+'_'+x+'.done/FrustrationData/LD')
		vectLD=np.zeros(int(resn[len(resn)-1])+1)	
		for ldline in locald.readlines():
			ldline=ldline.rstrip('\n')
			sp=ldline.split()
			if sp[0] != 'Res1':
				vectLD[int(sp[0])]=sp[1]
		
		locald.close()
		a=int(resn[0])
		vectmax=np.zeros(int(resn[len(resn)-1])+1)
		vecttotal=np.zeros(int(resn[len(resn)-1])+1)
		
		fst=open(jobID+'/Frustration/FrustrationResults/'+pdb+'_'+x+'.done/FrustrationData/'+pdb+'_'+x+'.pdb_mutational','r')	

		for line in fst.readlines():
			line=line.rstrip('\n')
			spfst=line.split()
			if spfst[0] != 'Res1' and sp[0] != 'NA' and sp[1] != 'NA':
				vecttotal[int(spfst[0])]+=1
				vecttotal[int(spfst[1])]+=1
				if spfst[13] == 'highly':	
					vectmax[int(spfst[0])]+=1
					vectmax[int(spfst[1])]+=1	
		fst.close()
		vectFV=np.zeros(int(resn[len(resn)-1])+1)
		for j in range(0,len(vectmax)):
			if vecttotal[j] !=0:
				vectFV[j] = (float(vectmax[j])/float(vecttotal[j]))
		outfst=open(jobID+'/Data/'+pdb+'_'+x+'.fstdata','w')
		for h in range(0,len(pomf)):
			outfst.write(str(resn[h])+' '+str(nocont[h])+' '+str(vectLD[h])+' '+str(vectFV[h])+' '+str(pohf[h])+' '+str(ponf[h])+' '+str(pomf[h])+'\n')
		outfst.close()
		
	return chains
	
def runmodel(JobId,PdbId,x):
	model = xgb.XGBClassifier()
	model.load_model('/home/maria/Documentos/FrustraPocket/savedModel.json')   
	columns = ['ID','NoC','PoHFC','LD','FV','PoNFC','PoMFC']
	protein_dataset = pd.read_csv(JobId+"/Data/"+PdbId+"_"+x+".fstdata", sep=" ", header=None, names=columns, index_col=False)
	a=protein_dataset['ID']
	protein_dataset = protein_dataset.drop(columns=['ID'])
	predictions = model.predict(protein_dataset)
	results = model.evals_result()
	results_dataset = protein_dataset.assign(ID=a,Class = predictions)
	results_dataset.to_csv(JobId+"/Data/"+PdbId+"_"+x+".csv")
	
	
def Pocket(JobId,PdbId,ch,tr):
	columns = ['W','NoC','PoHFC','LD','FV','PoNFC','PoMFC','ID','Class']
	data = pd.read_csv(JobId+"/Data/"+PdbId+"_"+ch+".csv", sep=",", names=columns, index_col=False)
	pockets=[]
	aux=[]
	p=0
	for x in range(1,len(data['Class'])):
		if int(data['Class'][x]) == 1:
			aux.append(data['ID'][x])
		else:
			if len(aux)>int(tr):
				p+=1
				a=aux.copy()
				pockets.append(a)
			aux.clear()
	return p,pockets
	
def output(JobId,PdbId,ch):
	p,pockets=Pocket(JobId,PdbId,ch,5)
	allcenter=open(JobId+'/OutPutFiles/CenterOfMass_'+ch,'w')
	if p == 0:
		p,pockets=Pocket(JobId,PdbId,ch,2)
	pdb=open(JobId+'/Frustration/'+PdbId+'_'+ch+'.pdb')
	lpdb=pdb.readlines()
	colores=['red', 'lightorange', 'lightblue ', 'paleyellow', 'sand','raspberry','purple', 'salmon', 'warmpink', 'yellow', 'hotpink','green', 'splitpea', 'palegreen', 'limon', 'blue','wheat', 'magenta',  'pink',  'aquamarine', 'orange', 'brightorange', 'deepolive', 'palecyan', 'lightpink']
	pdb.close()
	pymol=open(JobId+'/OutPutFiles/'+PdbId+'_'+ch+'.pml','w')
	pymol.write('load ../Frustration/'+PdbId+'_'+ch+'.pdb, protein\ncolor gray, protein\n')
	for i, values in enumerate(pockets):
		a=i+1
		out=open(JobId+'/OutPutFiles/Pocket_'+str(a)+'_'+ch+'.pdb','w')
		pymol.write('load Pocket_'+str(a)+'_'+ch+'.pdb, p'+str(a)+'\nshow surface,p'+str(a)+'\ncolor '+colores[i]+', p'+str(a)+'\n')
		for val in values:
			for res in lpdb:
				resa = res[22:26]
				for x in resa:
					if x.isalpha():
						resa.replace(x,'')
				for x in val:
					if x.isalpha():
						val.replace(x,'')
				if int(resa) == int(val):
					out.write(res)
		out.close()
	
		p=JobId+'/OutPutFiles/Pocket_'+str(a)+'_'+ch+'.pdb'
		calc = np.genfromtxt(p, skip_header=1, usecols=[6, 7, 8])
		v = calc.mean(axis=0)
		allcenter.write(p+' '+str(v[0])+' '+str(v[1])+' '+str(v[2])+'\n')
	allcenter.close()
	pymol.write('zoom all\n')
	pymol.close()	
