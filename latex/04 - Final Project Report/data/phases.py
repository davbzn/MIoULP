import csv
from math import *
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

plt.rc('text', usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('font',**{'family':'serif','serif':['Palatino']})
plt.rc('font', family='serif')

	####  data 	####
xfr0 = np.genfromtxt("fr0.csv", delimiter=',')
ypha = np.genfromtxt("phase.csv", delimiter=',')
ycor = np.genfromtxt("corr_phase.csv", delimiter=',')

# Creo un grafico la dimensione è in pollici
fig1 = plt.figure(figsize=(7, 4))
# Titolo del grafico
fig1.suptitle(r'Phase correction', y=0.98, fontsize=15)

######
# GRAFICO 1
f1 = fig1.add_subplot(1, 1, 1)
phase, = f1.plot(xfr0, ypha, '-', c='orange')
corrected, = f1.plot(xfr0, ycor, '-', c='red')

f1.set_ylabel(r'Phase', labelpad=2, fontsize=14)
f1.set_xlabel(r'Optical Frequency [PHz]', labelpad=2, fontsize=14)

#f1.set_ylim((-0.1, 1.1))
#f1.set_yticks(np.arange(-90, 1, 15))
#f1.set_xlim((-1.75,1.75))

f1.grid(True)

f1.legend([phase, corrected], ('Original\nphase', 'Corrected\n phase'), loc='lower left')  

fig1.subplots_adjust(left=0.12, right=0.96, top=0.92, bottom=0.12, hspace=0.085, wspace=0.05)
# mostra grafico
plt.show()
