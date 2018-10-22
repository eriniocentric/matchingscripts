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

for ifile in np.arange(1,10):

    file="3dhsttncutall_"+str(ifile)
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
group_by_name=mastertable.group_by('ID')
ascii.write(group_by_name,'comparetncut.matches',overwrite=True)

matches = ascii.read('compare3dhst.matches')

plt.axis([0,35,0,35])
plt.xlabel('SN2D')
plt.ylabel('SN1 (red) SN2 (blue)')
plt.plot([0,1,50],[0,1,50])
found = np.isin(mastertable['ID'], matches['ID'])
plt.plot(mastertable['SN2D'][found],mastertable['SN1'][found],color='black',marker='o',linestyle='None')
plt.plot(mastertable['SN2D'][found],mastertable['SN2'][found],color='black',marker='o',linestyle='None')

plt.plot(mastertable['SN2D'],mastertable['SN1'],'r.')
plt.plot(mastertable['SN2D'],mastertable['SN2'],'b.')

#plt.show()

plt.savefig('tncutsn_vs_sn2d.png')

plt.close()

plt.axis([0,35,0,35])
plt.xlabel('SN2D')
plt.ylabel('sqrt( SN1^2 + SN2^2 )')
plt.plot([0,1,50],[0,1,50])
found = np.isin(mastertable['ID'], matches['ID'])
snquad=np.sqrt(mastertable['SN1']**2+mastertable['SN2']**2)
plt.plot(mastertable['SN2D'][found],snquad[found],color='black',marker='o',linestyle='None')
plt.plot(mastertable['SN2D'],snquad,'b.')
plt.savefig('SNquad_vs_SNfit.png')
plt.close()


#plot histogram of objects from force photometry and those found as a detection in the source catalog

#open up forced photometry catalog                                                     
forced=ascii.read('all_output_gt8.cat')
plt.hist(forced['SN2D'],histtype='bar',alpha=0.5,color='orange',bins='auto',label='3dhst detections using 3dhst ra/dec')

plt.hist(matches['SN2D'],histtype='bar',alpha=0.5,color='blue',bins='auto',label='3dhst matches to HX detect catalog')
uniq=group_by_name.groups.indices

group_mean=group_by_name.groups.aggregate(np.mean)
plt.hist(group_mean['SN2D'],histtype='bar',alpha=0.5,color='green',bins='auto',label='3dhst matches with t*.ncut catalogs')
plt.xlabel('SN from rsp1a2')
plt.ylabel('number of LAEs recovered')
plt.legend()
plt.savefig('3dhstmatcheshist.png')
plt.close()

#plot up comparison of flux, ra, dec for HX catalog matches

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax1.scatter(matches['SN2D'],3600.0*(matches['RA']-matches['HXRA']))
ax1.set_xlabel('SN2D')
ax1.set_ylabel('RA-HXRA')
ax1.set_xlim(3,70)
ax1.set_ylim(-5,7)
ax2.scatter(matches['SN2D'],3600.0*(matches['DEC']-matches['HXDEC']))
ax2.set_xlabel('SN2D')
ax2.set_ylabel('DEC-HXDEC')
ax2.set_xlim(3,70)
ax2.set_ylim(-5,7)
ax3.scatter(matches['SN2D'],matches['HXFLUX']/matches['FLUX'])
ax3.set_xlabel('SN2D')
ax3.set_ylabel('HXFLUX/FLUX')
ax3.set_xlim(3,70)
ax3.set_ylim(-1,4)
ax4.scatter(matches['deltawave'],matches['separation'])
ax4.set_xlabel('deltawave (A) ')
ax4.set_ylabel('separation (arcsec)')
ax4.set_xlim(-1,3)
ax4.set_ylim(-1,6)

plt.savefig('3dhstscatterplots.png')
