#!/usr/bin/python

import sys
import os
import os.path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from astropy.io import fits
from astropy.table import Table, hstack, vstack, Column, MaskedColumn
from astropy.io import ascii

#this is from all the t*.ncut files from cosmos, gn and egs
#had to break the file up

table = Table(names=('ID','RA','DEC','WAVEOUT','FLUX','SIGMA','SN2D'),dtype=('S16','f8','f8','f8','f8','f8','f8'))

tableHX = Table(names=('HXWAVE','HXRA','HXDEC','SN1','SN2'),dtype=('f8','f8','f8','f8','f8')) 

tablesep = Table(names=('separation','deltawave'),dtype=('f8','f8'))

for ifile in np.arange(1,9):

    file="3dhsttcutall_"+str(ifile)
    print file
    data=ascii.read(file)
    print "Finished reading file: "+str(file)
    sep = 6./3600.
    wavesep = 6
    data3dhst=ascii.read('all_output_gt8.cat')

    output=np.array([])
    i=0
    while i < len(data3dhst):
        #print i
        RA = data3dhst['RA'][i]
        DEC = data3dhst['DEC'][i]
        wave = data3dhst['WAVEOUT'][i]
        deltawave=abs(data['col5'] - wave)
        distance = ((data['col19']-RA)**2 + (data['col20']-DEC)**2)**(1./2.)
#        print np.min(distance)
        if any ((distance < sep) & (deltawave<wavesep)):
            sel = np.where((distance<sep) & (deltawave<wavesep))
            sel2=np.argmin(deltawave[sel])
#        print data[sel][sel2]
#        print data3dhst[i]
            table.add_row(data3dhst['ID','RA','DEC','WAVEOUT','FLUX','SIGMA','SN2D'][i])
            tableHX.add_row(data['col19','col20','col5','col27','col28'][sel][sel2])
            tablesep.add_row((3600.*distance[sel][sel2],deltawave[sel][sel2]))
            print data3dhst['ID','RA','DEC','WAVEOUT','FLUX','SIGMA','SN2D'][i]
            print data['col19','col20','col5','col27','col28'][sel][sel2]
            print 3600.*distance[sel][sel2]
            print deltawave[sel][sel2]
         
        i += 1

mastertable=hstack([table,tableHX,tablesep])

plt.axis([0,35,0,35])
plt.xlabel('SN2D')
plt.ylabel('SN1 (red) SN2 (blue)')
plt.plot([0,1,50],[0,1,50])
plt.plot(mastertable['SN2D'],mastertable['SN1'],'r.')
plt.plot(mastertable['SN2D'],mastertable['SN2'],'b.')
plt.savefig('tncutsn_vs_sn2d.png')
