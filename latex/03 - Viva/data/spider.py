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
xfr0 = np.genfromtxt("spi_fr.csv", delimiter=',')
yspe = np.genfromtxt("spi_int.csv", delimiter=',')
yphi = np.genfromtxt("spi_phy.csv", delimiter=',')

# Creo un grafico la dimensione è in pollici
fig1 = plt.figure(figsize=(7, 4))
# Titolo del grafico
fig1.suptitle(r'Spectrum and Phase of $U_r$ from SEA-F-SPIDER', y=0.97, fontsize=15)

######
# GRAFICO 1
f1 = fig1.add_subplot(1, 1, 1)
spectrum, = f1.plot(xfr0, yspe/max(yspe), '-', c='blue')

f1.set_ylabel(r'Intensity [a.u.]', labelpad=2, fontsize=14)
#f1.set_xlabel(r'Optical Frequency [PHz]', labelpad=2, fontsize=14)
f1.set_xlabel(r'Optical Frequency [PHz]', labelpad=2, fontsize=14)

f1.set_ylim((-0.1, 1.27))
#f1.set_yticks(np.arange(-90, 1, 15))
f1.set_xlim((0.2,0.8))

f1.grid(True)

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

f2.set_xlabel('Wavelength [nm]')

f3 = f1.twinx()
phase, = f3.plot(xfr0, yphi, '-', c='orange')
f3.set_ylabel(u'Phase', labelpad=2, fontsize=14)

f1.legend([spectrum, phase],('Spectrum', 'Phase'), loc='center right')

#fig1.subplots_adjust(left=0.08, right=0.72, top=0.8, bottom=0.12, hspace=0.085, wspace=0.05)
fig1.subplots_adjust(left=0.08, right=0.94, top=0.8, bottom=0.12, hspace=0.085, wspace=0.05)
# mostra grafico
plt.show()
