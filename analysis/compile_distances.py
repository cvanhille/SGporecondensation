import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import os

if len(sys.argv) < 3:
	print("Error! Not enough arguments. Should be two: 1) File to analyse 2) Plot keyword (y/n, yes/no, true/false, True/False or plot/noplot)")
	exit()

file = sys.argv[1]
plot = sys.argv[2]
if plot == 'y' or plot == 'plot' or plot == 'yes' or plot == 'True' or plot == 'true':
	plot = True
elif plot == 'n' or plot == 'noplot' or plot == 'no' or plot == 'False' or plot == 'false':
	plot = False
else:
	print("Error! Second argument is for plotting. Should be yes or no (y/n, yes/no, true/false, True/False or plot/noplot)")
	exit()

d = pd.read_csv(file, sep = '\t', skiprows = 1, header = None, index_col = 0)

time = []
dist = []
dnor = []
dnum = []
drog = []
prog = []
pnor = []

for f in d.index:
	frame = d.loc[f,:]
	pore = np.array(frame[1].split(' ')).astype('float')
	porepos = pore[:3]
	porenor = pore[3:6]
	porer = pore[6]
	for i in frame.index[1:]:
		dr = frame[i]
		if isinstance(dr, str):
			drop = np.array(dr.split(' ')).astype(float)
			droppos = drop[1:4]
			dropn = drop[0]
			dropr = drop[4]
			distvec = droppos-porepos
			distvec[distvec >= 30] = distvec[distvec >= 30]-60
			distvec[distvec < -30] = distvec[distvec < -30]+60
			distnor = np.dot(distvec, porenor)
			time.append(f)
			dist.append(np.linalg.norm(distvec))
			dnor.append(distnor)
			dnum.append(dropn)
			drog.append(dropr)
			prog.append(porer)
			pnor.append(np.dot(porenor, np.array([1,0,0])))

time = np.array(time)
dist = np.array(dist)
dnor = np.array(dnor)
dnum = np.array(dnum)
drog = np.array(drog)
prog = np.array(prog)
pnor = np.array(pnor)

path = '%s/NR3.33/compiled_distances'%(file.split('NR3.33/')[0])
name = '%s/%s_%s_%s.txt'%(path, file.split('NR3.33/')[1].split('/')[0], file.split('NR3.33/')[1].split('/')[1], file.split('NR3.33/')[1].split('/')[2].split('_')[0])

if not os.access(path, os.F_OK):
	r = os.system('mkdir %s'%(path))

data = pd.DataFrame(index = np.arange(len(time)), columns = ['Time', 'Distance', 'Distance +', 'Droplet size', 'Droplet radius', 'Pore radius', 'Pore X normal'])
data['Time'] = time
data['Distance'] = dist
data['Distance +'] = dnor
data['Droplet size'] = dnum
data['Droplet radius'] = drog
data['Pore radius'] = prog
data['Pore X normal'] = pnor
data.to_csv(name)

if plot:
	print("Plotting...")

	plt.figure(figsize = (8,6))
	plt.plot(time, dist, 'o')
	plt.tick_params(axis = 'both', labelsize = 20)
	plt.xlabel('Time', fontsize = 25)
	plt.ylabel(r'Distance to pore $\Delta r$ [$\sigma$]', fontsize = 25)
	plt.tight_layout()
	plt.show()

	plt.figure(figsize = (8,6))
	plt.plot(time, dnor, 'o')
	plt.tick_params(axis = 'both', labelsize = 20)
	plt.xlabel('Time', fontsize = 25)
	plt.ylabel(r'Normal distance to pore $\Delta r_+$ [$\sigma$]', fontsize = 25)
	plt.tight_layout()
	plt.show()

	plt.figure(figsize = (8,6))
	plt.plot(time, dnum, 'o')
	plt.tick_params(axis = 'both', labelsize = 20)
	plt.xlabel('Time', fontsize = 25)
	plt.ylabel('Droplet size', fontsize = 25)
	plt.tight_layout()
	plt.show()

	plt.figure(figsize = (8,6))
	plt.plot(time, drog, 'o')
	plt.tick_params(axis = 'both', labelsize = 20)
	plt.xlabel('Time', fontsize = 25)
	plt.ylabel(r'Droplet radius $r_\mathrm{drop}$ [$\sigma$]', fontsize = 25)
	plt.tight_layout()
	plt.show()

	plt.figure(figsize = (8,6))
	plt.plot(time, dist/drog, 'o')
	plt.tick_params(axis = 'both', labelsize = 20)
	plt.xlabel('Time', fontsize = 25)
	plt.ylabel(r'$\Delta r / r_\mathrm{drop} $', fontsize = 25)
	plt.tight_layout()
	plt.show()

	plt.figure(figsize = (8,6))
	plt.plot(time, dnor/drog, 'o')
	plt.tick_params(axis = 'both', labelsize = 20)
	plt.xlabel('Time', fontsize = 25)
	plt.ylabel(r'$\Delta r_+ / r_\mathrm{drop} $', fontsize = 25)
	plt.tight_layout()
	plt.show()

	plt.figure(figsize = (8,6))
	plt.plot(time, prog, 'o')
	plt.tick_params(axis = 'both', labelsize = 20)
	plt.xlabel('Time', fontsize = 25)
	plt.ylabel(r'Pore radius $r_\mathrm{pore}$ [$\sigma$]', fontsize = 25)
	plt.tight_layout()
	plt.show()

	plt.figure(figsize = (8,6))
	plt.plot(time, drog/prog, 'o')
	plt.tick_params(axis = 'both', labelsize = 20)
	plt.xlabel('Time', fontsize = 25)
	plt.ylabel(r'$r_\mathrm{drop} / r_\mathrm{pore}$', fontsize = 25)
	plt.tight_layout()
	plt.show()

	plt.figure(figsize = (8,6))
	plt.plot(time, pnor, 'o')
	plt.tick_params(axis = 'both', labelsize = 20)
	plt.xlabel('Time', fontsize = 25)
	plt.ylabel(r'Pore normal $X$ component', fontsize = 25)
	plt.tight_layout()
	plt.show()
