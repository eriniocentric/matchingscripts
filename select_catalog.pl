#accesses the SDSS MPA-JHU DR7 catalog http://www.mpa-garching.mpg.de/SDSS/DR7/
#select galaxies from the SDSS that starforming, Halpha bright and are in the SDP H-Atlas Field at 9h30m,+0d,30arcmin

use PDL;
use PDL::NiceSlice;

#open info Catalog
$galINFO = rfits 'Research/DustTemp/gal_info_dr7_v5_2.fit';

$ra = $galINFO->{'RA'};#in degrees
$dec = $galINFO->{'DEC'};#in degrees
$z = $galINFO->{'Z'}; #redshift
$zflag = $galINFO->{'Z_WARNING'};# if non-zero, redshift is not worthy

@spectrotype = $galINFO->{'SPECTROTYPE'};#array ref containing spectrotype of 'GALAXY', 'AGN'

#open total SFR catalog
$grandSFR = rfits 'Research/DustTemp/gal_totsfr_dr7_v5_2.fits';

$logsfr = $grandSFR->{'AVG'};#use the average SFR from the PDF calculated in Brinchmann+2004


$sel = which ( #($ra < 15 | $ra > 354) & #$ra > 129 & $ra < 141 &
	       $zflag < 1 &
	       $z < 0.08 &
	       $logsfr > log10(50) & 
	       $logsfr < 3.9);

wcols $ra($sel),$dec($sel),'Research/DustTemp/coords.dat';  

wcols $ra($sel),$dec($sel),$z($sel),$logsfr($sel),'Research/DustTemp/selected_cat.dat';
