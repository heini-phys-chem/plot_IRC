# plot_IRC
Python programs to extract IRC energies from gaussien output file and plots them. <br />
The first program (get_IRC_data.py) takes forward and reverse IRC paths and orders them correctly. The output are 2 cPickle binaries containing energies and xyz's. Then plot_irc.py will load both cPickle files (or as many different IRC's as you want) and plots them. <br />

needed libraries:
- matplotlib
- cclib
- numpy
- cPickle

command line should look like:
```
python3 get_IRC_data.py 'file'.log
python plot_irc.py
```
