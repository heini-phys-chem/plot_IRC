#!/usr/bin/env python3

import sys, os
import numpy as np
import cclib

import matplotlib.pyplot as plt
import seaborn as sns


sns.set_style('whitegrid')
conv = 23. # Ha to eV

NAME = {1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 10:'Ne', 11:'Na', 12:'Mg', 13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl', 35:'Br', 53:'I'}


def read_data(f):
	data = cclib.io.ccread(f)

	return data

def extract_data(data):

	coords   = data.atomcoords
	atomnums = data.atomnos
	labels   = [NAME[atomnum] for atomnum in atomnums]
	numAtoms = len(labels)
	energies = data.scfenergies

	return numAtoms, labels, coords, energies

def check_pos(coords):
	if np.linalg.norm(coords[0] - coords[3]) < np.linalg.norm(coords[1] - coords[2]):
		return False 
	else:
		return True

def write_xyz(f, coords, labels, numAtoms, energy):

	file = open(f, 'a')
	file.write("{}\n{}\n".format(numAtoms, energy))

	for i, atom in enumerate(coords):
		file.write("{}\t{}\t{}\t{}\n".format(labels[i], atom[0], atom[1], atom[2]))

	file.close()

if __name__ == "__main__":

	#os.system("rm -f *.xyz")

	filename = sys.argv[1]

	data                               = read_data(filename)
	numAtoms, labels, coords, energies = extract_data(data)

	e_ts   = np.max(energies)
	idx_ts = np.argmax(energies)
	e1     = np.min(energies)
	idx_e1 = np.argmin(energies)

	for i in range(1, len(energies)):
		if np.abs(energies[i-1]-energies[i]) > .85:
			e2     = energies[i-1]
			idx_e2 = i-1

	f_ts    = filename[:-4] + "_ts.xyz"
	f_react = filename[:-4] + "_react.xyz"
	f_prod  = filename[:-4] + "_prod.xyz"

	write_xyz(f_ts, coords[idx_ts], labels, numAtoms, energies[idx_ts])

	is_prod = check_pos(coords[idx_e1])

	if is_prod == True:
		if conv*(np.abs(energies[idx_e1] - energies[idx_e2])) != 0.0:
			write_xyz(f_prod, coords[idx_e1], labels, numAtoms, energies[idx_e1])
			write_xyz(f_react, coords[idx_e2], labels, numAtoms, energies[idx_e2])
			print(f_ts, conv*(e_ts - energies[idx_e2]), conv*(np.abs(energies[idx_e1] - energies[idx_e2])))
	else:
		if conv*(np.abs(energies[idx_e1] - energies[idx_e2])) != 0.0:
			write_xyz(f_prod, coords[idx_e2], labels, numAtoms, energies[idx_e2])
			write_xyz(f_react, coords[idx_e1], labels, numAtoms, energies[idx_e1])
			print(f_ts, conv*(e_ts - energies[idx_e1]), conv*(np.abs(energies[idx_e1] - energies[idx_e2])))
