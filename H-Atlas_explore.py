
import numpy as np
import astropy
import aplpy
import matplotlib.pyplot as plt

from astropy.io import fits

hatlasfile = fits.open('Research/H-ATLAS/DR1/HATLAS_DR1_CATALOGUE_V1.2.FITS')
hdata = hatlasfile[1].data

#route data into easy to deal with variables
ID = hdata.field('IDNAME')
RA = hdata.field('RA')
DEC = hdata.field('DEC')

zpec = hdata.field('Z_SPEC')
zqual = hdata.field('Z_QUAL') #only spectroscopic redshifts with Z_QUAL>=3 are recommended to be reliable
sdssidhatlas = hdata.field('SDSS_OBJID') #Optical IDs from SDSS DR7

#create photometry data array

wave = [3.368,4.618,12.082,22.194,100,160,250,350,500]
phot = np.array([hdata.field('W1_FLUX_GAMA'),hdata.field('W2_FLUX_GAMA'),hdata.field('W3_FLUX_GAMA'),hdata.field('W4_FLUX_GAMA'),hdata.field('F100BEST'),hdata.field('F160BEST'),hdata.field('F250BEST'),hdata.field('F350BEST'),hdata.field('F500BEST')])


#open SDSS DR7 catalogs to get SDSS IDs and logSFRs

sdssidfile = fits.open('Research/DustTemp/CAS-IDs.fits')
sdssid = sdssidfile[1].data.field('OBJID')
sdssidfile.close()

sdssfile = fits.open('Research/DustTemp/gal_info_dr7_v5_2.fit')
sdssdata = sdssfile[1].data
sdssfile.close()

sdssRA = sdssdata.field('RA')
sdssDEC = sdssdata.field('DEC')
sdsszpec = sdssdata.field('Z')
sdsszflag = sdssdata.field('Z_WARNING')
sdssgaltype = sdssdata.field('SPECTROTYPE') #galaxy type either QSO or GALAXY

sdssfile2 = fits.open('Research/DustTemp/gal_totsfr_dr7_v5_2.fits')
sdssdata2 = sdssfile2[1].data
sdssfile2.close()

sfr = sdssdata2.field('AVG')

#find sdss sfr to print to new catalog

sfr_hatlas = list()

for x in np.nditer(sdssidhatlas):
    indx = np.where(x == sdssid)[0]
    sfr_hatlas.append(sfr[indx])

#update new H-ATLAS Fits table with SFR info

newhdu = fits.BinTableHDU.from_columns([fits.Column(name='SFR', format = 'E', array=sfr)])

hatlasfile.append(newhdu)
hatlasfile.writeto('Research/H-ATLAS/DR1/HATLAS_DR1_addedSFR.fits')
