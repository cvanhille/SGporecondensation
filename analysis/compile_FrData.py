import glob
import os
import sys

if len(sys.argv) < 2:
	print("Error! Not enough arguments. Should be two: 1) seed for glob.glob from .../NR3.33)")
	exit()

seed = sys.argv[1]

files = glob.glob('/nfs/scistore15/saricgrp/cvanhill/sims/SG_pore_condensation/vesicle/NR3.33/%s'%(seed))
gpath = '/nfs/scistore15/saricgrp/cvanhill/sims/SG_pore_condensation/vesicle/NR3.33/compiled_FrData'
if not os.access(gpath, os.F_OK):
	r = os.system('mkdir %s'%(gpath))

for file in files:
	sfile = file.split('NR3.33/')[1]
	name = '%s_%s_%s.txt'%(sfile.split('/')[0], sfile.split('/')[1].split('rac_yes')[0], sfile.split('yes/FrData/')[1].split('_')[0])
	r = os.system('cp %s %s/%s'%(file, gpath, name))