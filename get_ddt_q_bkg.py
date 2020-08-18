# - - - - - - - - - - 
# Code to convert 3rd part of IDL zeropointing code to python3
# - - - - - - - - 
import numpy as np
import astropy.units as u
import scipy.io as scio
import matplotlib.pyplot as plt
import astropy.coordinates as coord

from astroquery.vizier import Vizier
from astropy.table import Table,Column
from make_jh_sav import cat_match
from itertools import zip_longest
Vizier.ROW_LIMIT = 99999
# - - - - - - - - 
# Input Files:
# - SerpensSouth_ALL-PHOT.csv
# - SerpensSouth_W-PHOT
# - serpenssouth_sexcat_w.fits

# Output Files:
# - SerpensSouth_WJH.csv
# - - - - - - - - 

# - - - - -  F U N C T I O N S - - - - -
#
# FUNCTION to calculate zeropoint
#
def find_zeropoint(phottab,wtab,wsexcat):

	sexcat_match, w_match, sexcat_ind, w_ind = cat_match(wsexcat,wtab,1.0*u.arcsec,'ALPHA_J2000','DELTA_J2000','rd','dd',CONVERT=False)
	emko = 1.85
	qmmko = -0.043

	# Determine zeropoint by assuming q of background stars
	zpwmko0 = (phottab['jm'][w_ind]+((1+emko)*2.5*np.log10(wsexcat[sexcat_ind]['FLUX_AUTO']))+(emko*phottab['hm'][w_ind])-qmmko)/(1+emko)
	goodzpmko = np.where((np.isfinite(zpwmko0))&(wsexcat[sexcat_ind]['FLAGS'] == 0)&(np.isfinite(phottab['jme'][w_ind]))&(np.isfinite(phottab['hme'][w_ind]))&(wsexcat[sexcat_ind]['CLASS_STAR'] > 0.5))	
	# Iterative fit
	for i in range(5):
		zpmko = np.nanmedian(zpwmko0[goodzpmko])
		zpmkoerr = np.nanmedian(abs(zpwmko0[goodzpmko]-zpmko))
		print(zpmko,zpmkoerr,len(goodzpmko[0]))
		goodzpmko = np.where((np.isfinite(zpwmko0))&(wsexcat[sexcat_ind]['FLAGS'] == 0)&(np.isfinite(phottab['jme'][w_ind]))&(np.isfinite(phottab['hme'][w_ind]))&(wsexcat[sexcat_ind]['CLASS_STAR'] > 0.5)&(wsexcat[sexcat_ind]['FLUXERR_AUTO']/wsexcat[sexcat_ind]['FLUX_AUTO'] < 0.1)&(abs(zpwmko0-zpmko) < 4.5*zpmkoerr))
 
	# Determine zeropoint from simulated W band mags using all available photmetry
	zpwsyn0 = (wtab['wcest'][w_ind] + 2.5*np.log10(wsexcat[sexcat_ind]['FLUX_AUTO']))
	goodzpsyn = np.where((np.isfinite(zpwsyn0))&(wtab['wceste'][w_ind] < 0.108)&(wsexcat[sexcat_ind]['FLAGS'] == 0)&(wsexcat[sexcat_ind]['CLASS_STAR'] > 0.5)&(wsexcat[sexcat_ind]['FLUXERR_AUTO']/wsexcat[sexcat_ind]['FLUX_AUTO'] < 0.1))
	# Iterative fit
	for i in range(5):	
		zpsyn = np.median(zpwsyn0[goodzpsyn])
		zpsynerr = np.median(abs(zpwsyn0[goodzpsyn]-zpsyn))
		print(zpsyn,zpsynerr,len(goodzpsyn[0]))
	goodzpsyn = np.where((np.isfinite(zpwsyn0))&(wtab['wceste'][w_ind] < 0.108)&(wsexcat[sexcat_ind]['FLAGS'] == 0)&(wsexcat[sexcat_ind]['CLASS_STAR'] > 0.5)&(wsexcat[sexcat_ind]['FLUXERR_AUTO']/wsexcat[sexcat_ind]['FLUX_AUTO'] < 0.1)&(abs(zpwsyn0-zpsyn) < 4.5*zpsynerr))

	wm = zpsyn-2.5*np.log10(wsexcat['FLUX_AUTO'])
	wme = 1.085*wsexcat['FLUXERR_AUTO']/wsexcat['FLUX_AUTO']

	# Write to file
	# Columns: RA,Dec,W,We,J,Je,H,He,WS,JS,HS
	final_tab = Table()
	ra = Column(np.around(np.array(phottab['rd'][w_ind]),6),name='RA')
	dec = Column(np.around(np.array(phottab['dd'][w_ind]),6),name='Dec')
	w = Column(np.around(np.array(wm[sexcat_ind]),2),name='W')
	we = Column(np.around(np.array(wme[sexcat_ind]),2),name='WERR')
	j = Column(np.around(np.array(phottab['jm'][w_ind]),2),name='J')
	je = Column(np.around(np.array(phottab['jme'][w_ind]),2),name='JERR')
	h = Column(np.around(np.array(phottab['hm'][w_ind]),2),name='H')
	he = Column(np.around(np.array(phottab['hme'][w_ind]),2),name='HERR')
	final_tab.add_columns([ra,dec,w,we,j,je,h,he])
	final_tab.write('SerpensSouth_WJHpy.dat',format='ascii',overwrite=True)

	return

# - - - - - - - - 
def main():

	# Read in photometry tables
	phottab = Table.read('SerpensSouth_ALL-PHOT.csv',format='csv')
	wtab = Table.read('SerpensSouth_W-PHOT.csv',format='csv')
	wsexcat = Table.read('serpenssouth_sexcat_w.fits',format='fits')

	find_zeropoint(phottab,wtab,wsexcat) 

	return

main()
