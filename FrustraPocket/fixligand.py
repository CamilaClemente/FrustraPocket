import os
import sys
import os.path as path

pdb=open(sys.argv[1]+'/'+sys.argv[2]+'_aux.pdbqt','r')
out=open(sys.argv[1]+'/'+sys.argv[2]+'.pdbqt','w')

for l in pdb.readlines():
	l=l[:-1]
	if len(l) > 18:
		if l[18] == ' ':
			out.write(l[0]+l[1]+l[2]+l[3]+l[4]+l[5]+l[6]+l[7]+l[8]+l[9]+l[10]+l[11]+l[12]+l[13]+l[14]+l[15]+l[16]+l[17]+l[19]+l[20]+l[21]+l[22]+l[23]+l[24]+l[25]+l[26]+l[27]+l[28]+l[29]+l[30]+l[31]+l[32]+l[33]+l[34]+l[35]+l[36]+l[37]+l[38]+l[39]+l[40]+l[41]+l[42]+l[43]+l[44]+l[45]+l[46]+l[47]+l[48]+l[49]+l[50]+l[51]+l[52]+l[53]+l[54]+l[55]+l[56]+l[57]+l[58]+l[59]+l[60]+l[61]+l[62]+l[63]+l[64]+l[65]+l[66]+l[67]+l[68]+l[69]+l[70]+l[71]+l[72]+l[73]+l[74]+l[75]+l[76]+l[77]+l[78]+'\n')
		else:
			out.write(l+'\n')
	else:
			out.write(l+'\n')
pdb.close()
out.close()

