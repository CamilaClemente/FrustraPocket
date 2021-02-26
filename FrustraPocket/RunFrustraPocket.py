import os
import sys
import os.path as path
import numpy as np

def Chains(direc, pdb):
	chains=[]
	#https://files.rcsb.org/header/1NRE.pdb	
	wget='wget \'https://files.rcsb.org/header/'+pdb+'.pdb\' -O '+direc+'/'+pdb+'header.pdb'
	os.system(wget)
	grep='grep \"CHAIN: \" '+direc+'/'+pdb+'header.pdb > '+direc+'/aux'
	os.system(grep)
	rm='rm '+direc+'/'+pdb+'header.pdb'
	os.system(rm)
	aux=open(direc+'/aux','r')
	for line in aux.readlines():
		sp=line.split()
		for i in range(3,len(sp)):
			c=sp[i]
			if c[-1] == ";" or c[-1] == ",":
				c=c[:-1]
			t=splitPDB(direc, pdb, c)
			if t == 1:
				chains.append(c)
	aux.close()
	rm='rm '+direc+'/aux'
	os.system(rm)
	return chains

def splitPDB(direc, pdb, chain):
	pdbo=open(direc+'/'+pdb+'.pdb','r')
	c=0
	r=0
	for lpdb in pdbo.readlines():
		atom=lpdb[0]+lpdb[1]+lpdb[2]+lpdb[3]
		if atom == 'ATOM' and lpdb[21] == chain:
			c+=1
	if c > 150:
		r=1
	pdbo.close()
	return r
	
def Fstandlden(dfrustra,dchain,pdb,chain):
	vect=[]
	frust=open(dfrustra+'FrustrationData/'+pdb+'_'+chain+'.pdb_mutational_5adens','r')
	out=open(dchain+'/'+pdb+'_aux.fstdata','w')
	d=1
	for lineade in frust.readlines():
		spadens=lineade.split()
		if d==1:
			d=2
		else:
			out.write(spadens[0]+' '+spadens[6]+' ')
			auxld=open(dfrustra+'FrustrationData/'+pdb+'_'+chain+'.pdb_mutational','r')
			for lld in auxld.readlines():
				spld=lld.split()
				if spadens[0] == spld[0]:
					out.write(spld[4])
					break
				if spadens[0] == spld[1]:
					out.write(spld[5])
					break
			out.write('\n')
			auxld.close()
	frust.close()
	out.close()
	sort='cd '+dchain+'/;sort -u -r -k2 '+pdb+'_aux.fstdata > '+pdb+'.fstdata;rm '+pdb+'_aux.fstdata'
	os.system(sort)
	

def FrustaPocket (fit,ldt,dchain,dfrustra,pdb,chain): # frustration index threshold (fit) and local density threshold (ldt), pockets directory results, frustraresultdirectory and pdb id
	vect=[]# empty vector for residues 
	out=open(dchain+'/'+pdb+'.pockets','w')
	frust=open(dchain+'/'+pdb+'.fstdata','r')
	pocket=1
	for linefrus in frust.readlines():
		d=0
		spfrust=linefrus.split()
		if spfrust[0] == 'NA':
			print('NA')
		else:
			if len(spfrust) == 3:
				if float(spfrust[1]) >= fit and float(spfrust[2]) >= ldt:
					if len(vect) > 0:
						d=0
						for i in range (0, len(vect)):
							if vect[i] == spfrust[0]:
								d=1
					if d == 0: 
						auxmut=open(dfrustra+'FrustrationData/'+pdb+'_'+chain+'.pdb_mutational','r')
						vectaux=[]
						vectaux.append(spfrust[0])
						for mutline in auxmut.readlines():
							spmut=mutline.split()
							if spmut[0] == spfrust[0]:
								vectaux.append(spmut[1])
							elif spmut[1] == spfrust[0]: 
								vectaux.append(spmut[0])
						auxmut.close()
						if pocket == 1 and len(vectaux) > 1:
							out.write('Pocket '+str(pocket))
							for res in range(0,len(vectaux)):
								out.write(' '+vectaux[res])
							out.write('\n')
							pocket+=1
							vect=vect+vectaux
						elif len(vectaux) > 1:
							vaux=[]
							vi=0
							for i in range(0,len(vectaux)):
								r=0
								for j in range(0,len(vect)):
									if vectaux[i] == vect[j]:
										r=1
								if r==0:
									vaux.append(vectaux[i])
							if len(vectaux) == len(vaux):
								out.write('Pocket '+str(pocket))
								for res in range(0,len(vectaux)):
									out.write(' '+vectaux[res])
								out.write('\n')
								pocket+=1
								vect=vect+vectaux
	out.close()
	frust.close()
	return pocket
	

#----- Creating Directories ----

pdb=sys.argv[1]

pipedir=os.getcwd()+'/'
direc=pipedir+'job.'+pdb
fst=direc+'/'+pdb+'.done/FrustrationData/'
mkdir='mkdir '+direc
rm='rm -r '+direc
os.system(rm)
os.system(mkdir)
mkdir='mkdir '+direc+'/Pockets/'
os.system(mkdir)
#----- Download PDB file -----

pathPDB=os.getcwd()+'/'+pdb+'.pdb'

if path.exists(pathPDB):
	cp='cp '+pathPDB+' '+direc+'/'+pdb+'_aux.pdb'
	os.system(cp)
else:
	wget='wget \'http://www.rcsb.org/pdb/files/'+pdb+'.pdb\' -O '+direc+'/'+pdb+'_aux.pdb'
	os.system(wget)

fixpdb='python3 fixpdb.py '+direc+' '+pdb
os.system(fixpdb)
#----- FrustrAR -----
chains = Chains(direc, pdb)
dchain=''

for i in range(0,len(chains)):
	dchain=direc+'/'+chains[i]
	mkdir='mkdir '+dchain
	os.system(mkdir)
	outfst=open(dchain+'/frustrar.R','w')
	outfst.write('library(frustratometeR)\nOrderList <- c("'+pdb+'.pdb")\nPdb_mut <- dir_frustration(PdbsDir = \"'+direc+'/\",Chain = \"'+chains[i]+'\", OrderList = OrderList, Mode = "mutational", ResultsDir = \"'+dchain+'/\")')
	outfst.close()
	rsc='cd '+dchain+';Rscript frustrar.R'
	os.system(rsc)
	#----- Generating Pockets -----
	dfrustra=dchain+'/'+pdb+'_'+chains[i]+'.done/'
	Fstandlden(dfrustra,dchain,pdb,chains[i])
	pocket=FrustaPocket (0.13,2.6,dchain,dfrustra,pdb,chains[i])
	
	if int(pocket) <= 3:
		pocket=FrustaPocket (0.2,1,dchain,dfrustra,pdb,chains[i])

outpml=open(direc+'/Pockets/VisualizationPockets.pml','w')
outpml.write('load '+pdb+'.pdb\nshow cartoon, all\n')

cp='cp '+direc+'/'+pdb+'.pdb '+direc+'/Pockets/'+pdb+'.pdb'
os.system(cp)
allcenter=open(direc+'/Pockets/Pocket_centerofmass.txt','w')

for j in range(0,len(chains)):
	pockets=open(direc+'/'+chains[j]+'/'+pdb+'.pockets','r')
	for lpock in pockets.readlines():
		lpock=lpock[:-1]
		sppock=lpock.split()
		outp=open(direc+'/Pockets/'+sppock[1]+'_'+chains[j]+'.pdb','w')
		mayor=0
		for l in range(2, len(sppock)):
			if int(sppock[l]) < 0:
				sppock[l]=int(sppock[l])*-1
			grep='grep \"'+str(sppock[l])+' '+chains[j]+'\" '+direc+'/'+chains[j]+'/'+pdb+'_'+chains[j]+'.done/FrustrationData/'+pdb+'_'+chains[j]+'.pdb_mutational_5adens > '+direc+'/a'
			os.system(grep)
			aa=open(direc+'/a','r')
			laa=aa.readline()
			spa=laa.split()
			if float(mayor) < float(spa[6]):
				mayor=spa[6]
			aa.close()
		#p=score/(len(sppock)-2)
		outp.write('REMARK	'+str(mayor)+'\n')
		
		rm='rm '+direc+'/a'
		os.system(rm)
		outpml.write('load '+sppock[1]+'_'+chains[j]+'.pdb\nshow surface,'+sppock[1]+'_'+chains[j]+'\n')
		for l in range(2, len(sppock)):
			pdbn=open(direc+'/'+chains[j]+'/'+pdb+'_'+chains[j]+'.done/FrustrationData/'+pdb+'_'+chains[j]+'.pdb','r')
			for lpdb in pdbn.readlines():
				if len(lpdb) > 53:
					aa=lpdb[22]+lpdb[23]+lpdb[24]+lpdb[25]
					aa=int(aa)
					if aa == int(sppock[l]):
						outp.write(lpdb)
			pdbn.close()
		outp.close()
		p=direc+'/Pockets/'+sppock[1]+'_'+chains[j]+'.pdb'
		#p='\"'+p+'\"'
		#print(p)
		a = np.genfromtxt(p, skip_header=1, usecols=[6, 7, 8])
		v = a.mean(axis=0)
		allcenter.write(p+' '+str(v[0])+' '+str(v[1])+' '+str(v[2])+'\n')
allcenter.close()
outpml.write('zoom all')
outpml.close()

if len(sys.argv) > 2:
	awk='cd '+direc+'/;awk \'{if ($1 == "ATOM"){print} if ($1=="TER"){print}}\' '+pdb+'.pdb > '+pdb+'_clean.pdb'
	os.system(awk)
	lig=sys.argv[2]
	pligand='cp '+pipedir+'prepare_ligand4.py '+direc+'/prepare_ligand4.py'
	os.system(pligand)
	precep='cp '+pipedir+'prepare_receptor4.py '+direc+'/prepare_receptor4.py'
	os.system(precep)
	cplig='cp '+pipedir+lig+'.pdb '+direc+'/'+lig+'.pdb'
	os.system(cplig)

	pythonpath='/home/XXX/MGLTools-1.5.7rc1/bin/pythonsh'
	preceptor='cd '+direc+'/;'+pythonpath+' prepare_receptor4.py -A hydrogens -r '+pdb+'_clean.pdb -o '+pdb+'_aux.pdbqt'
	fixpdbqt='python3 fixpdbqt.py '+direc+' '+pdb
	pliga='cd '+direc+'/;'+pythonpath+' prepare_ligand4.py -A hydrogens -l '+lig+'.pdb -o '+lig+'_aux.pdbqt'
	fixligand='python3 fixligand.py '+direc+' '+lig
	os.system(pliga)
	os.system(preceptor)
	os.system(fixpdbqt)
	os.system(fixligand)
	rm='cd '+direc+'/;rm '+pdb+'_aux.pdbqt '+lig+'_aux.pdbqt'
	os.system(rm)
	center=open(direc+'/Pockets/Pocket_centerofmass.txt','r')
	outvina=open(direc+'/VinaResults.txt','w')
	for lines in center.readlines():
		lines=lines[:-1]
		splines=lines.split()
		print('Vina Results for pocket:'+splines[0])
		vina='cd '+direc+'/;vina --receptor '+pdb+'.pdbqt --ligand '+lig+'.pdbqt --center_x  '+splines[1]+' --center_y '+splines[2]+' --center_z '+splines[3]+' --size_x 30 --size_y 30 --size_z 30 --out '+splines[0]+'_ligand'
		print(vina)
		os.system(vina)
		vhead='cd '+direc+'/;head -2 '+splines[0]+'_ligand > aux'
		os.system(vhead)
		vaux=open(direc+'/aux','r')
		vline=vaux.readlines()
		outvina.write(splines[0]+' '+vline[1])
	center.close()
	outvina.close()
