import sys
import numpy as np
import cPickle
import matplotlib.pyplot as plt
#import cclib
#from cclib.parser import ccopen
import seaborn as sns
sns.set(color_codes=True)

def get_irc(coords):
  y = []
  for i,mol in enumerate(coords):
    y.append(np.linalg.norm(mol[0]-mol[3]) - np.linalg.norm(mol[1]-mol[5]))

  return(y)
  #plt.scatter(y,energies)
  #plt.plot(y, energies, color='black')
  #plt.show()

if __name__ == "__main__":

  calc_en = cPickle.load(open('irc_reference_dft_energies.cpickle', 'rb'))
  calc_mol = cPickle.load(open('irc_reference_dft_mols.cpickle', 'rb'))
  ml_en = cPickle.load(open('irc_lqa_2_energies.cpickle', 'rb'))
  ml_mol = cPickle.load(open('irc_lqa_2_mols.cpickle', 'rb'))
  ml2_en = cPickle.load(open('irc_lqa_4_energies.cpickle', 'rb'))
  ml2_mol = cPickle.load(open('irc_lqa_4_mols.cpickle', 'rb'))

  calc_irc = get_irc(calc_mol)
  ml_irc = get_irc(ml_mol)
  ml2_irc = get_irc(ml2_mol)

  # plot
  fig = plt.figure
  ax = plt.subplot(111)
  # scatter
  ax.scatter(calc_irc,calc_en, s=35, label="QM")
  ax.scatter(ml_irc,ml_en, color='r', s=35, label="LQA (step size 2)")
  ax.scatter(ml2_irc,ml2_en, color='g', s=35, label="LQA (step size 4)")
  ax.legend(fontsize=25)

  # lines
  #ax.plot(calc_irc, calc_en, color='black')
  #ax.plot(ml_irc, ml_en, color='black')
  #ax.plot(ml2_irc, ml2_en, color='black')

  ax.tick_params(labelsize=25)

  ax.set_xlabel("IRC", fontsize=25)
  ax.set_ylabel("Energy [kcal/mol]", fontsize=25)

  plt.show()

