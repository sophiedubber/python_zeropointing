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
def find_zeropoint(phottab,wphottab,w):

	wphottab = wphottab.filled(np.nan)
	w_match, wphottab_match, ww_ind, syn_ind = cat_match(w,wphottab,1.0*u.arcsec,'ALPHA_J2000','DELTA_J2000','rd','dd',CONVERT=False)
	emko = 1.85
	qmmko = -0.043

	# Determine zeropoint by assuming q of background stars
	zpwmko0 = (phottab['jm'][syn_ind]+((1+emko)*2.5*np.log10(w[ww_ind]['FLUX_AUTO']))+(emko*phottab['hm'][syn_ind])-qmmko)/(1+emko)
	goodzpmko = np.where((np.isfinite(zpwmko0))&(w[ww_ind]['FLAGS'] == 0)&(np.isfinite(phottab['jme'][syn_ind]))&(np.isfinite(phottab['hme'][syn_ind]))&(w[ww_ind]['CLASS_STAR'] > 0.5))	
	# Iterative fit
	for i in range(5):
		zpmko = np.nanmedian(zpwmko0[goodzpmko])
		zpmkoerr = np.nanmedian(abs(zpwmko0[goodzpmko]-zpmko))
		print(zpmko,zpmkoerr,len(goodzpmko[0]))
		goodzpmko = np.where((np.isfinite(zpwmko0))&(w[ww_ind]['FLAGS'] == 0)&(np.isfinite(phottab['jme'][syn_ind]))&(np.isfinite(phottab['hme'][syn_ind]))&(w[ww_ind]['CLASS_STAR'] > 0.5)&(w[ww_ind]['FLUXERR_AUTO']/w[ww_ind]['FLUX_AUTO'] < 0.1)&(abs(zpwmko0-zpmko) < 4.5*zpmkoerr))
 
	# Determine zeropoint from simulated W band mags using all available photmetry
	zpwsyn0 = (wphottab['wcest'][syn_ind] + 2.5*np.log10(w[ww_ind]['FLUX_AUTO']))
	goodzpsyn = np.where((np.isfinite(zpwsyn0))&(wphottab['wceste'][syn_ind] < 0.108)&(w[ww_ind]['FLAGS'] == 0)&(w[ww_ind]['CLASS_STAR'] > 0.5)&(w[ww_ind]['FLUXERR_AUTO']/w[ww_ind]['FLUX_AUTO'] < 0.1))
	# Iterative fit
	for i in range(5):	
		zpsyn = np.median(zpwsyn0[goodzpsyn])
		zpsynerr = np.median(abs(zpwsyn0[goodzpsyn]-zpsyn))
		print(zpsyn,zpsynerr,len(goodzpsyn[0]))
		goodzpsyn = np.where((np.isfinite(zpwsyn0))&(wphottab['wceste'][syn_ind] < 0.108)&(w[ww_ind]['FLAGS'] == 0)&(w[ww_ind]['CLASS_STAR'] > 0.5)&(w[ww_ind]['FLUXERR_AUTO']/w[ww_ind]['FLUX_AUTO'] < 0.1)&(abs(zpwsyn0-zpsyn) < 4.5*zpsynerr))

	flagtest = np.where(w['FLAGS'][ww_ind] == 0)
	wm = zpsyn-2.5*np.log10(w['FLUX_AUTO'])
	wme = 1.085*w['FLUXERR_AUTO']/w['FLUX_AUTO']

	# Write to file
	# Columns: RA,Dec,W,We,J,Je,H,He,WS,JS,HS
	final_tab = Table()
	ra = Column(np.around(np.array(w['ALPHA_J2000'][ww_ind[flagtest]]),6),name='RA')
	dec = Column(np.around(np.array(w['DELTA_J2000'][ww_ind[flagtest]]),6),name='Dec')
	w = Column(np.around(np.array(wm[ww_ind[flagtest]]),2),name='W')
	we = Column(np.around(np.array(wme[ww_ind[flagtest]]),2),name='WERR')
	j = Column(np.around(np.array(phottab['jm'][syn_ind[flagtest]]),2),name='J')
	je = Column(np.around(np.array(phottab['jme'][syn_ind[flagtest]]),2),name='JERR')
	h = Column(np.around(np.array(phottab['hm'][syn_ind[flagtest]]),2),name='H')
	he = Column(np.around(np.array(phottab['hme'][syn_ind[flagtest]]),2),name='HERR')
	av = Column(np.around(np.array(wphottab['avest'][syn_ind[flagtest]]),2),name='AV')
	final_tab.add_columns([ra,dec,w,we,j,je,h,he,av])
	#final_tab.write('SerpensSouth_WJHpy.dat',format='ascii',overwrite=True)

	return final_tab

# - - - - - - - - 
def main():

	# Read in photometry tables
	phottab = Table.read('SerpensSouth_ALL-PHOT.csv',format='csv')
	wphottab = Table.read('SerpensSouth_W-PHOT_edit.csv',format='csv')
	w = Table.read('serpenssouth_sexcat_w.fits',format='fits')

	find_zeropoint(phottab,wphottab,w) 

	return 

main()
