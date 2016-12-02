from pylab import *
from numpy import *
from sys import * 
import csv

''' 
    This script extracts the IRC Path of a Gaussian output file using matplotlib
'''

# variable list
x       = []
y       = []
data    = []

# check if input file is specified
if len(argv) < 2:
    print '\n##############################'
    print 'command line should look like:'
    print 'python plot.py <file>.log'
    print '###############################\n'
    sys.exit()

script, file_in = argv

# open file, read in data set
with open(file_in) as f:
    copy = False
    for line in f:
        # start tag
        if line.strip() == 'Reaction path calculation complete.':
            copy = True
        # end tag
        elif line.strip() == 'IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC':
            copy = False
        # store data in list without blamk lines    
        elif copy:
            line = filter(lambda x: not x == '\n', line)
            data.append(line)

# total number of calculated points on Rx
points = int(data[-3][-3:]) + 1

# extrax x and y values from data set
for val in data[6:6+points]:
    y.append(float(val[23:30]))
    x.append(float(val[32:40]))

# reverse list and convert from Hartree to kcal/mol (*627.509)
y = [i*627.509 for i in list(reversed(y))]

# react, prod and TS energy - after reordering rel energies and after unit convertion
react   = y[0]
prod    = y[-1]
TS      = data[1][-11:]

# generate a figure
fig = figure()
fig.suptitle('Minimal energy Path for Claisen Rxn', fontsize=14, fontweight='bold')
ax = fig.add_subplot(111)
fig.subplots_adjust(top=.9)

# set subtitle and label axes
ax.set_title('TS energy = ' + TS + ' Hartree')
ax.set_xlabel('Reaction coordinate')
ax.set_ylabel('Energy rel. to TS [kcal/mol]')

# text input
text_in = ('''$E_{react}=$\t %r kcal/mol 
$E_{prod}=$\t %r kcal/mol''' % (react,prod))

# text in plot
ax.text(10,-5, text_in, fontsize=14,  bbox={'facecolor':'none', 'alpha':0.5, 'pad':10})

# generate plot and set axes
plot(x,y,'o',x,y,'k')
xlim([x[0]-2,x[-1]+2])
ylim([min(y)-2,max(y)+2])

# show plot
show()
