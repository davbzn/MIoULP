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
xtau = np.genfromtxt("tau.csv", delimiter=',')
yint = np.genfromtxt("intensity.csv", delimiter=',')

# Creo un grafico la dimensione è in pollici
fig1 = plt.figure(figsize=(7, 4))
# Titolo del grafico
fig1.suptitle(r"Intensity of pulse", y=0.97, fontsize=15)

######
# GRAFICO 1
f1 = fig1.add_subplot(1, 1, 1)
#f2.set_xscale('log')

interference = f1.errorbar(x=xtau, y=yint, fmt='.:', c='blue')

#fase100 = f1.errorbar(x=f100, y=ph100, fmt='.:', c='green')

#teo = f1.errorbar(x=np.logspace(10,6E5,500), y=180/pi*np.arctan(np.logspace(10,6E5,500)/(1+200000*11)/8), fmt='.:', c='green')

f1.set_ylabel(r'Intensity [a.u.]', labelpad=2, fontsize=14)
f1.set_xlabel(r'Delay $\tau$ [fs]', labelpad=2, fontsize=14)

#f1.text(2E4, -37.5, 'G=101x', #r'$\nu_0$', size=12, va='center', ha='center')
#f1.text(2E5, -37.5, 'G=11x', size=12, va='center', ha='center')

f1.set_ylim((0, 1.1))
#f1.set_yticks(np.arange(-90, 1, 15))
#f1.set_xlim((5,6E5))

#plt.setp(f1.get_xticklabels(), visible=False)

f1.grid(True)
#f1.legend((fase10, fase100), ("Gain = 11x", "Gain = 101x"), 'lower left', prop={'size': 12})  
######

# questo imposta i bordi del grafico
fig1.subplots_adjust(left=0.09, right=0.97, top=0.91, bottom=0.13, hspace=0.085, wspace=0.05)
# mostra grafico
plt.show()
