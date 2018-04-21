import sys
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
import cclib

from cclib.parser import ccopen

NAME = {
    1:  'H'  ,
    2:  'He' ,
    3:  'Li' ,
    4:  'Be' ,
    5:  'B'  ,
    6:  'C'  ,
    7:  'N'  ,
    8:  'O'  ,
    9:  'F'  ,
   10:  'Ne' ,
   11:  'Na' ,
   12:  'Mg' ,
   13:  'Al' ,
   14:  'Si' ,
   15:  'P'  ,
   16:  'S'  ,
   17:  'Cl' ,
   18:  'Ar' ,
   19:  'K'  ,
   20:  'Ca' ,
   21:  'Sc' ,
   22:  'Ti' ,
   23:  'V'  ,
   24:  'Cr' ,
   25:  'Mn' ,
   26:  'Fe' ,
   27:  'Co' ,
   28:  'Ni' ,
   29:  'Cu' ,
   30:  'Zn' ,
   31:  'Ga' ,
   32:  'Ge' ,
   33:  'As' ,
   34:  'Se' ,
   35:  'Br' ,
   36:  'Kr' ,
   37:  'Rb' ,
   38:  'Sr' ,
   39:  'Y'  ,
   40:  'Zr' ,
   41:  'Nb' ,
   42:  'Mo' ,
   43:  'Tc' ,
   44:  'Ru' ,
   45:  'Rh' ,
   46:  'Pd' ,
   47:  'Ag' ,
   48:  'Cd' ,
   49:  'In' ,
   50:  'Sn' ,
   51:  'Sb' ,
   52:  'Te' ,
   53:  'I'  ,
   54:  'Xe' ,
   55:  'Cs' ,
   56:  'Ba' ,
   57:  'La' ,
   58:  'Ce' ,
   59:  'Pr' ,
   60:  'Nd' ,
   61:  'Pm' ,
   62:  'Sm' ,
   63:  'Eu' ,
   64:  'Gd' ,
   65:  'Tb' ,
   66:  'Dy' ,
   67:  'Ho' ,
   68:  'Er' ,
   69:  'Tm' ,
   70:  'Yb' ,
   71:  'Lu' ,
   72:  'Hf' ,
   73:  'Ta' ,
   74:  'W'  ,
   75:  'Re' ,
   76:  'Os' ,
   77:  'Ir' ,
   78:  'Pt' ,
   79:  'Au' ,
   80:  'Hg' ,
   81:  'Tl' ,
   82:  'Pb' ,
  82:  'Pb' ,
  83:  'Bi' ,
  84:  'Po' ,
  85:  'At' ,
  86:  'Rn' ,
  87:  'Fr' ,
  88:  'Ra' ,
  89:  'Ac' ,
  90:  'Th' ,
  91:  'Pa' ,
  92:  'U'  ,
  93:  'Np' ,
  94:  'Pu' ,
  95:  'Am' ,
  96:  'Cm' ,
  97:  'Bk' ,
  98:  'Cf' ,
  99:  'Es' ,
 100:  'Fm' ,
 101:  'Md' ,
 102:  'No' ,
 103:  'Lr' ,
 104:  'Rf' ,
 105:  'Db' ,
 106:  'Sg' ,
 107:  'Bh' ,
 108:  'Hs' ,
 109:  'Mt' ,
 110:  'Ds' ,
 111:  'Rg' ,
 112:  'Cn' ,
 114:  'Uuq',
 116:  'Uuh'}

def merge(l1, l2):
  yield from l1
  yield from l2

def split_list(lst, mols, mid):
  """ split IRC list in half and rearrange to correct order
      also leave away first geometry (input geometry)
  """
  firstHalf = lst[mid:]
  secondHalf = lst[:mid]
  mols1 = mols[mid:]
  mols2 = mols[:mid]

  if firstHalf[-1] > secondHalf[-1]:
    new_list = list(merge(reversed(firstHalf), secondHalf))
    new_mols = list(merge(reversed(mols1), mols2))
  else:
    new_list = list(merge(reversed(secondHalf), firstHalf))
    new_mols = list(merge(reversed(mols2), mols1))

  return new_list, new_mols

def plot_irc(energies,coords):
  y = []
  for i,mol in enumerate(coords):
    y.append(np.linalg.norm(mol[0]-mol[3]) - np.linalg.norm(mol[1]-mol[5]))

  plt.scatter(y,energies)
  plt.plot(y, energies, color='black')
  plt.show()

if __name__ == "__main__":
  # read in log file
  filename = sys.argv[1]
  base = filename[:-4]

  mylogfile = ccopen(filename)
  data = mylogfile.parse()

  # get number of atoms, coordinates, energies, atom numbers and transition state energy
  nAtoms = data.natom
  mols = data.atomcoords
  ts = mols[0]
  mols = mols[:-1]
  energies = data.scfenergies
  atoms = [NAME[i] for i in data.atomnos]
  ts_energy = max(energies)

  y = range(len(energies))
  plt.scatter(y,energies)
  plt.plot(y, energies, color='black')
  plt.show()


  plot_irc(energies,mols)

  diff = []
  for i,energy in enumerate(energies):
    if i == 0: continue
    energy_before = energies[i-1]
    if (energy*627.509)-(energy_before*627.509) > 100:
      mid=i

  to_plot, mols_new = split_list(energies, mols, mid)

  new_list = np.asarray(to_plot)
  new_list -= new_list[0]
  new_list *= 23

  plot_irc(new_list, mols_new)

  print("\n save data")
  with open(str(base) + "_energies.cpickle", 'wb') as f:
    cPickle.dump(new_list, f, protocol=2)
  with open(str(base) + "_mols.cpickle", 'wb') as f:
    cPickle.dump(mols_new, f, protocol=2)

  for mol in mols_new:
    print(nAtoms)
    print("")
    for i,atom in enumerate(mol):
      print(str(atoms[i]), *atom)
  #print(nAtoms)
  #print("")
  #for i,atom in enumerate(ts):
  #  print(str(atoms[i]), *atom)

