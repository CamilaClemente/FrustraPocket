import os
import sys

pdb=open(sys.argv[1]+'/'+sys.argv[2]+'_aux.pdb','r')
out=open(sys.argv[1]+'/'+sys.argv[2]+'.pdb','w')

for lpdb in pdb.readlines():
	atom=lpdb[0]+lpdb[1]+lpdb[2]+lpdb[3]
	if atom == 'ATOM':
		if lpdb[26] == ' ':
			out.write(lpdb)
		else:
			out.write(lpdb[0]+lpdb[1]+lpdb[2]+lpdb[3]+lpdb[4]+lpdb[5]+lpdb[6]+lpdb[7]+lpdb[8]+lpdb[9]+lpdb[10]+lpdb[11]+lpdb[12]+lpdb[13]+lpdb[14]+lpdb[15]+lpdb[16]+lpdb[17]+lpdb[18]+lpdb[19]+lpdb[20]+lpdb[21]+lpdb[22]+lpdb[23]+lpdb[24]+lpdb[25]+' '+lpdb[27]+lpdb[28]+lpdb[29]+lpdb[30]+lpdb[31]+lpdb[32]+lpdb[33]+lpdb[34]+lpdb[35]+lpdb[36]+lpdb[37]+lpdb[38]+lpdb[39]+lpdb[40]+lpdb[41]+lpdb[42]+lpdb[43]+lpdb[44]+lpdb[45]+lpdb[46]+lpdb[47]+lpdb[48]+lpdb[49]+lpdb[50]+lpdb[51]+lpdb[52]+lpdb[53]+lpdb[54]+lpdb[55]+lpdb[56]+lpdb[57]+lpdb[58]+lpdb[59]+lpdb[60]+lpdb[61]+lpdb[62]+lpdb[63]+lpdb[64]+lpdb[65]+lpdb[66]+lpdb[67]+lpdb[68]+lpdb[69]+lpdb[70]+lpdb[71]+lpdb[72]+lpdb[73]+lpdb[74]+lpdb[75]+lpdb[76]+lpdb[77]+'\n')
	else:
		out.write(lpdb)
pdb.close()
out.close()	
rm='rm '+sys.argv[1]+'/'+sys.argv[2]+'_aux.pdb'	
os.system(rm)
