import csv
from math import *
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import math

plt.rc('text', usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
plt.rc('font', family='serif')

def getWaveLength(frequency):
    return 299.792458/frequency

	####  data 	####
xfr0 = np.genfromtxt("fr0.csv", delimiter=',')
xwle = np.genfromtxt("wlen.csv", delimiter=',')
yspe = np.genfromtxt("spectrum.csv", delimiter=',')
ygau = np.genfromtxt("gauss.csv", delimiter=',')

# Creo un grafico la dimensione è in pollici
fig1 = plt.figure(figsize=(7, 4))
# Titolo del grafico
fig1.suptitle(r'Spectrum and Gaussian filter', y=0.98, fontsize=15)

######
# GRAFICO 1
f1 = fig1.add_subplot(1, 1, 1)
filtered, = f1.plot(xfr0, ygau*yspe/max(yspe[550:800]), '-', c='magenta')
spectrum, = f1.plot(xfr0, yspe/max(yspe[550:800]), '-', c='blue')
gauss,    = f1.plot(xfr0, ygau, '-', c='red')

f1.set_ylabel(r'Intensity [a.u.]', labelpad=2, fontsize=14)
#f1.set_xlabel(r'Optical Frequency [PHz]', labelpad=2, fontsize=14)
f1.set_xlabel(r'Optical Frequency [PHz], Wavelength [nm]', labelpad=2, fontsize=14)

f1.set_ylim((-0.1, 1.1))
#f1.set_yticks(np.arange(-90, 1, 15))
f1.set_xlim((0.1,0.7))

f1.grid(True)

f1.legend([spectrum, gauss, filtered],('Spectrum', 'Gauss\n filter', 'Spectrum\n filtered'),
	  loc='center right', bbox_to_anchor=(1.375, 0.5) )#, prop={'size': 12})  
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

f2 = f1.twiny()  # ax2 is responsible for "top" axis and "right" axis

# get the primary axis x tick locations in plot units
xtickloc = f1.get_xticks() 
# set the second axis ticks to the same locations
f2.set_xticks(xtickloc)
# calculate new values for the second axis tick labels, format them, and set them
x2labels = [r'{:.4g}'.format(math.ceil(x/10)*10) for x in getWaveLength(xtickloc)]
f2.set_xticklabels(x2labels)
# force the bounds to be the same
f2.set_xlim(f1.get_xlim()) 

f1.text(0.76, 1.156, r'[nm]', size=12, va='center', ha='center')

#f2.set_xlabel('Wavelength [nm]')

#fig1.subplots_adjust(left=0.08, right=0.72, top=0.8, bottom=0.12, hspace=0.085, wspace=0.05)
fig1.subplots_adjust(left=0.08, right=0.75, top=0.86, bottom=0.12, hspace=0.085, wspace=0.05)
# mostra grafico
plt.show()