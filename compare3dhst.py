#!/usr/bin/python

import sys
import os
import os.path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn, hstack
from astropy.io import ascii

data=ascii.read('all3dhst.cat')

sep = 5./3600.
wavesep = 5
data3dhst=ascii.read('all_output_gt8.cat')

table = Table(names=('ID','RA','DEC','WAVEOUT','FLUX','SIGMA','SN2D'),dtype=('S16','f8','f8','f8','f8','f8','f8'))
tableHX = Table(names=('HXID','HXRA','HXDEC','HXWAVE','HXFLUX','HXSIGMA','HXSN2D'),dtype=('S16','f8','f8','f8','f8','f8','f8')) 
tablesep = Table(names=('separation','deltawave'),dtype=('f8','f8'))

i=0
while i < len(data3dhst):
    RA = data3dhst['RA'][i]
    DEC = data3dhst['DEC'][i]
    wave = data3dhst['WAVEOUT'][i]
    deltawave=abs(data['WAVE'] - wave)
    distance = ((data['RA']-RA)**2 + (data['DEC']-DEC)**2)**(1./2.)
#    print np.min(distance)
    if any ((distance < sep) & (deltawave<wavesep)):
        sel = np.where((distance<sep) & (deltawave<wavesep))
        sel2=np.argmin(deltawave[sel])
#        print data[sel][sel2]
#        print data3dhst[i]
        table.add_row(data3dhst['ID','RA','DEC','WAVEOUT','FLUX','SIGMA','SN2D'][i])
        tableHX.add_row(data['ID','RA','DEC','WAVE','FLUX','SIGMA','SN2D'][sel][sel2])
        tablesep.add_row((3600.*distance[sel][sel2],deltawave[sel][sel2]))
        
        print data3dhst['ID','RA','DEC','WAVEOUT','FLUX','SIGMA','SN2D'][i]
        print data['ID','RA','DEC','WAVE','FLUX','SIGMA','SN2D'][sel][sel2]
         
    i += 1

mastertable=hstack([table,tableHX,tablesep])
ascii.write(mastertable,'compare3dhst.matches',overwrite=True)




