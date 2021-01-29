import numpy as np
import numpy as np
import sys

def unPBC(x, boxLx,flag):
	if (flag=='ortho'):
		if (x > boxLx*0.5):
			x = x - boxLx
		elif (x <= -boxLx*0.5):
			x = x + boxLx
		return x

def initialize_coords_idnum_qdat():
	global coords,idnum,bonds,angles,totatoms,gro_dat,qdat,torsions
	coords = {}
	idnum = {}
	bonds = {}
	angles = {}
	qdat = {}
	for i in list_atoms:
		coords[i]= np.zeros((0,3))
		idnum[i] = np.zeros((0,1),dtype=int)
		bonds[i] = np.zeros((0,2),dtype=int)
		angles[i] = np.zeros((0,3),dtype=int)
	torsions = np.zeros((0,4),dtype=int)
	totatoms=len(gro_dat)
	print "Total atoms = ",totatoms
	i=-1
	for line in gro_dat:
		i = i+1
		linesplit = line.split()
		ID = linesplit[1]
		x = float(linesplit[3])
		y = float(linesplit[4])
		z = float(linesplit[5])
		aid = int(linesplit[2])
		coords[ID] = np.append(coords[ID],[[x,y,z]],axis=0)
		idnum[ID] = np.append(idnum[ID],[[aid]],axis=0)
		qdat[aid]=q[ID]


def diff_mat_CN_calc():
	global diff_mat,bw0_id_mat,bw1_id_mat,diffwrtbw0,diffwrtbw1,wrtbw0,wrtbw1
	diff_mat = np.zeros((coords[bond_bw[0]].shape[0],coords[bond_bw[1]].shape[0]))
	bw0_id_mat = np.zeros((idnum[bond_bw[0]].shape[0],idnum[bond_bw[1]].shape[0]),dtype=int)
	bw1_id_mat = np.zeros((idnum[bond_bw[0]].shape[0],idnum[bond_bw[1]].shape[0]),dtype=int)
	row=-1
	for bond1 in coords[bond_bw[0]]:
		row=row+1
		col=-1
		for bond2 in coords[bond_bw[1]]:
			col=col+1
			dr = bond1-bond2
			dr[0] = unPBC(dr[0],xB,xf)
			dr[1] = unPBC(dr[1],yB,yf)
			dr[2] = unPBC(dr[2],zB,zf)
			dr2 = np.dot(dr,dr)
			diff_mat[row,col] = dr2
			bw0_id_mat[row,col] = int(idnum[bond_bw[1]][col])
			bw1_id_mat[row,col] = int(idnum[bond_bw[0]][row])
	
	print 'Shape of matrix %s (along row) vs %s (along cols) with id of %s' % (bond_bw[0],bond_bw[1],bond_bw[1]), bw0_id_mat.shape
	print 'Shape of matrix %s (along row) vs %s (along cols) with id of %s' % (bond_bw[0],bond_bw[1],bond_bw[0]), bw1_id_mat.shape
	print idnum[bond_bw[0]].shape,idnum[bond_bw[1]].shape
	print "Only row element should change"
	print bw0_id_mat
	print "Only col element should change"
	print  bw1_id_mat
	
	wrtbw0ind = diff_mat.argsort(axis=1)[0:diff_mat.shape[0],0:CN_bond[bond_bw[0]]]
	wrtbw1ind = diff_mat.argsort(axis=0)[0:CN_bond[bond_bw[1]],0:diff_mat.shape[1]]
	
	wrtbw0 = np.zeros(wrtbw0ind.shape,dtype=int)
	wrtbw1 = np.zeros(wrtbw1ind.shape,dtype=int)
	diffwrtbw0 = np.zeros(wrtbw0ind.shape)
	diffwrtbw1 = np.zeros(wrtbw1ind.shape)
	
	#print idnum[bond_bw[1]]
	
	row = 0
	while row < wrtbw0ind.shape[0]:
		wrtbw0[row,:] = bw0_id_mat[row,wrtbw0ind[row,:]]
		diffwrtbw0[row,:] = diff_mat[row,wrtbw0ind[row,:]]
		row = row+1
	col = 0
	while col < wrtbw1ind.shape[1]:
		wrtbw1[:,col] = bw1_id_mat[wrtbw1ind[:,col],col]
		diffwrtbw1[:,col] = diff_mat[wrtbw1ind[:,col],col]
		col = col+1

	#print diffwrtbw0.shape[0],wrtbw0.shape[0],idnum[bond_bw[0]].shape[0]

def check_diff_mat_CN_calc():
	global diff_mat,bw0_id_mat,bw1_id_mat,diffwrtbw0,diffwrtbw1,wrtbw0,wrtbw1,find_skip_index
	print np.sqrt(np.amin(diffwrtbw0)),np.sqrt(np.average(diffwrtbw0)),np.sqrt(np.amax(diffwrtbw0))
	print 'Checking any pairs that violate bonding distances wrt %s' % bond_bw[0]
	countatoms = 0
	checkOatoms1 = []
	for i,j,k in zip(diffwrtbw0,wrtbw0,idnum[bond_bw[0]]):
		for l,m in zip(i,j):
			#if (l > 0.0324):
			#	#print np.sqrt(l),k,m
			if (l < 0.0324):
				countatoms+=1
				checkOatoms1+=[m]
	checkOatoms1 = np.unique(checkOatoms1)
	print "Number of zeolite framework pairs wrt %s = %d" % (bond_bw[0],countatoms)
	print '\n'
	nloops=0
	print diffwrtbw1.shape[1],wrtbw1.shape[1],idnum[bond_bw[1]].shape[0]
	print np.sqrt(np.amin(diffwrtbw1)),np.sqrt(np.average(diffwrtbw1)),np.sqrt(np.amax(diffwrtbw1))
	print '\nChecking any pairs that violate bonding distances wrt %s' % bond_bw[1]
	find_skip_index = []
	countatoms=0
	checkOatoms2 = []
	for i,j,k in zip(np.transpose(diffwrtbw1),np.transpose(wrtbw1),idnum[bond_bw[1]]):
		flag = -1 #Checking for any O that is not bonded to any Si
		for l,m in zip(i,j):
			nloops = nloops+1
			if (l < 0.0324):
				flag = 0
				countatoms+=1
				checkOatoms2+=[k]
			#if (l > 0.0324):
			#	#print np.sqrt(l),k,m
		if (flag == -1):
			print "\n .. Found",k
			find_skip_index+=[k]
			print "\n .. Noted",k
	find_skip_index = np.array(find_skip_index)
	checkOatoms2 = np.unique(checkOatoms2)
	
	if (len(checkOatoms1)==len(checkOatoms2) and np.array_equiv(np.sort(checkOatoms1),np.sort(checkOatoms2))):
		print "The number of unique O atoms wrt Si and O are the same"
	else:
		print "The number of unique O atoms wrt Si and O are NOT same"
	print "Number of zeolite framework pairs wrt %s = %d" % (bond_bw[1],countatoms)
	
	print 'Checked %d in loop over all %s atoms <= This is it' % (nloops,bond_bw[1])

def old_to_new_gro():
	global gro_dat_new,gro_dat,find_skip_index
	print find_skip_index
	gro_dat_new = np.chararray(shape=(len(gro_dat)-len(find_skip_index)),itemsize=64)
	i = -1
	for line in gro_dat:
		linesplit = line.split()
		if (int(linesplit[2]) not in find_skip_index):
			i = i+1
			linesplit[2]=str(i+1)
			gro_dat_new[i] = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (int(linesplit[0][:-3]),linesplit[0][-3:],linesplit[1],int(linesplit[2]),float(linesplit[3]),float(linesplit[4]),float(linesplit[5]))
		else:
			print linesplit[2]
	gro_dat = np.chararray(shape=(len(gro_dat_new)),itemsize=64)
	gro_dat[:] = gro_dat_new[:]
	#print gro_dat

def add_H_bonds_angles():
	global idnum,coords,diffwrtbw0,diffwrtbw1,wrtbw0,wrtbw1,bonds,angles,torsions,CN_bond,bond_bw,xB,xf,yB,yf,zB,zf,gro_dat_new,totatoms,qdat
	print '\n Adding H to any %s with is uncoordinated' % bond_bw[1]
	for i in range(idnum[bond_bw[1]].shape[0]):
		k = idnum[bond_bw[1]][i]
		c = coords[bond_bw[1]][i]
		for j in range(diffwrtbw1[:,i].shape[0]):
			l = diffwrtbw1[j,i]
			m = wrtbw1[j,i]
			if (l > 0.0324 and k not in find_skip_index):
				rnew = coords[bond_bw[1]][i]
				direction = int(sys.argv[2])
				if c[direction]<0.5*box[direction]:
					rnew[direction] = c[direction]-0.0945 #0.1
				elif c[direction]>=0.5*box[direction]:
					rnew[direction] = c[direction]+0.0945 #0.1
				newID = totatoms+1
				wrtbw1[j,i] = newID
				coords['H'] = np.append(coords['H'],[rnew],axis=0)
				idnum['H'] = np.append(idnum['H'],[[newID]],axis=0)
				bonds[bond_bw[1]] = np.append(bonds[bond_bw[1]],[[int(k),int(newID)]],axis=0)
				angles[bond_bw[1]] = np.append(angles[bond_bw[1]],[[int(wrtbw1[0,i]),int(k),int(wrtbw1[1,i])]],axis=0)
				SIaid = wrtbw1[j-2+1,i]
				#print wrtbw1[:,i],SIaid,wrtbw1[j,i]
				SIindexidnum = np.where(idnum[bond_bw[0]]==SIaid)
				#print SIindexidnum,SIaid,idnum[bond_bw[0]][SIindexidnum]
				for OSI in wrtbw0[SIindexidnum[0],:][0]:
					#print SIaid,OSI,k,wrtbw0[SIindexidnum[0],:]
					if OSI != k:
						torsions = np.append(torsions, [[int(OSI),int(SIaid),int(k),int(wrtbw1[j,i])]],axis=0)
				dr = c-rnew
				dr[0] = unPBC(dr[0],xB,xf)
				dr[1] = unPBC(dr[1],yB,yf)
				dr[2] = unPBC(dr[2],zB,zf)
				dr2 = np.dot(dr,dr)
				diffwrtbw1[j,i] = dr2
				gro_dat_new = np.append(gro_dat_new,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (2,'ZSH','H',newID,rnew[0],rnew[1],rnew[2]))
				qdat[newID]=q['H']
				qdat[int(k)]=q['OH']
				totatoms=newID

def main():
	global list_atoms,q,mass,qdat,coords,idnum,bonds,angles,torsions,k_torsion,bond_eq,k_bond_eq,angle_eq,k_angle_eq,bond_bw,CN_bond,angle_bw,gro_dat,totatoms,box,xf,xB,yf,yB,zf,zB,diff_mat,bw0_id_mat,bw1_id_mat,wrtbw0,wrtbw1,diffwrtbw0,diffwrtbw1,find_skip_index
	list_atoms = ['Si','O','H']
	q = {}
	q['Si'] = '1.5'#'2.10000'
	q['O'] = '-0.75'#'-1.05000'
	q['H'] = '0.435'#'0.42500'
	q['OH'] = '-0.810'#'-0.95'
	mass = {}
	mass['Si'] = '28.08550'
	mass['O'] = '15.99940'
	mass['H'] = '1.00794'
	bond_eq = 0.0945
	k_bond_eq = 463700.0
	angle_eq = 108.5 #109.47
	k_angle_eq = 460.621 # kB=0.00831446261815324 kJ/mol/K #251.0 
	bond_bw = ['Si','O']
	CN_bond = {}
	CN_bond[bond_bw[0]]=4
	CN_bond[bond_bw[1]]=2
	angle_bw = ['Si','O','H']
	k_torsion = ['1.359914','4.079741','0','-5.43965','0','0']
	f = open(sys.argv[1])
	gro_dat = f.read().split('\n')
	f.close()
	totatoms=np.array(gro_dat[1]).astype(int)
	box = np.array(gro_dat[-2].split()).astype(float)
	print 'The box dimensions are [x,y,z] : ',box
	gro_dat = gro_dat[2:-2]
	i = -1
	xf = 'ortho'
	xB = box[0]
	yf = 'ortho'
	yB = box[1]
	zf = 'ortho'
	zB = box[2]

main()
iter = 0
while True:
	iter+=1
	print "This is iteration %d" % iter
	initialize_coords_idnum_qdat()
	diff_mat_CN_calc()
	check_diff_mat_CN_calc()
	old_to_new_gro()
	initialize_coords_idnum_qdat()
	diff_mat_CN_calc()
	check_diff_mat_CN_calc()
	if (len(find_skip_index)==0):
		break

add_H_bonds_angles()
			

#### Write atoms ####
top_dat = np.chararray(shape=(len(gro_dat_new),8),itemsize=8)
i = -1
for line in gro_dat_new:
	i = i+1
	linesplit = line.split()
	top_dat[i,0] = linesplit[2]
	top_dat[i,1] = linesplit[1]
	top_dat[i,2] = str(1)
	top_dat[i,3] = linesplit[0][1:]
	top_dat[i,4] = linesplit[1]
	top_dat[i,5] = linesplit[2]
	top_dat[i,6] = qdat[int(linesplit[2])]
	top_dat[i,7] = mass[linesplit[1]]
g = open(sys.argv[1].split('/')[-1].split('.')[0]+'_wH_atoms.top','w')
for r in range(np.size(top_dat,axis=0)):
	for c in range(np.size(top_dat,axis=1)):
		g.write(str('\t'))
		g.write(str(top_dat[r,c]))
	g.write('\n')

g.close()

### Write bonds ####
line = np.chararray(shape=(bonds[bond_bw[1]].shape[0],5),itemsize=9)
i=-1
for bond in bonds[bond_bw[1]]:
	i = i+1
	line[i,0] = bond[0]
	line[i,1] = bond[1]
	line[i,2] = str(1)
	line[i,3] = str(bond_eq)
	line[i,4] = ';'+str(k_bond_eq)

g = open(sys.argv[1].split('.')[0]+'_wH_bonds.top','w')
for r in range(np.size(line,axis=0)):
	for c in range(np.size(line,axis=1)):
		g.write(str('\t'))
		g.write(str(line[r,c]))
	g.write('\n')
g.close()

### Write angles ####
line = np.chararray(shape=(angles[bond_bw[1]].shape[0],6),itemsize=9)
i=-1
for angle in angles[bond_bw[1]]:
	i = i+1
	line[i,0] = angle[0]
	line[i,1] = angle[1]
	line[i,2] = angle[2]
	line[i,3] = str(1)
	line[i,4] = str(angle_eq)
	line[i,5] = str(k_angle_eq)

g = open(sys.argv[1].split('.')[0]+'_wH_angles.top','w')
for r in range(np.size(line,axis=0)):
	for c in range(np.size(line,axis=1)):
		g.write(str('\t'))
		g.write(str(line[r,c]))
	g.write('\n')
g.close()

### Write torsions ###
line = np.chararray(shape=(torsions.shape[0],11),itemsize=9)
i=-1
for torsionquad in torsions:
	i = i+1
	line[i,0] = torsionquad[0]
	line[i,1] = torsionquad[1]
	line[i,2] = torsionquad[2]
	line[i,3] = torsionquad[3]
	line[i,4] = '3'
	line[i,5] = k_torsion[0]
	line[i,6] = k_torsion[1]
	line[i,7] = k_torsion[2]
	line[i,8] = k_torsion[3]
	line[i,9] = k_torsion[4]
	line[i,10] = k_torsion[5]

g = open(sys.argv[1].split('.')[0]+'_wH_torsion.top','w')
for r in range(np.size(line,axis=0)):
	for c in range(np.size(line,axis=1)):
		g.write(str('\t'))
		g.write(str(line[r,c]))
	g.write('\n')
g.close()


#### Write gro ####
g = open(sys.argv[1].split('.')[0]+'_wH.gro','w')
g.write('Added H to %s \n' % sys.argv[1])
g.write('%d \n' % totatoms)
for line in gro_dat_new:
	g.write(line)
	g.write("\n")
for line in box:
	g.write(str(line))
	g.write("  ")
g.close()

"""
#### Write Bonds ####
line = np.chararray(shape=(wrtbw0.shape[0]*4,5),itemsize=9)
print len(idnum[bond_bw[0]]),wrtbw0.shape
i = -1
for id,neighbors in zip(idnum[bond_bw[0]],wrtbw0):
	for neighbor in neighbors:
		i = i+1
		line[i,0] = id[0]
		line[i,1] = neighbor
		line[i,2] = str(1)
		line[i,3] = str(bond_eq)
		line[i,4] = str(k_bond_eq)

g = open(sys.argv[1].split('.')[0]+'_bonds.top','w')
for r in range(np.size(line,axis=0)):
	for c in range(np.size(line,axis=1)):
		g.write(str('\t'))
		g.write(str(line[r,c]))
	g.write('\n')
g.close()


#### Write Angles ####
line = np.chararray(shape=(wrtbw0.shape[0]*6,6),itemsize=9)
print len(idnum[bond_bw[0]]),wrtbw0.shape
i = -1
for id,neighbors in zip(idnum[bond_bw[0]],wrtbw0):
	for n1 in range(len(neighbors)):
		for n2 in range(len(neighbors[n1+1:])):
			i = i+1
			line[i,0] = neighbors[n1]
			line[i,1] = id[0]
			line[i,2] = neighbors[n1+1+n2]
			line[i,3] = str(2)
			line[i,4] = str(angle_eq)
			line[i,5] = str(k_angle_eq)

g = open(sys.argv[1].split('.')[0]+'_angles.top','w')
for r in range(np.size(line,axis=0)):
	for c in range(np.size(line,axis=1)):
		g.write(str('\t'))
		g.write(str(line[r,c]))
	g.write('\n')
g.close()
"""
