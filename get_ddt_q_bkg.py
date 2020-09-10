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

	# Convert numerical spectral types to strings
	spt_list = np.around(np.array(wphottab['spt']),2)
	spt_str = []
	for i in range(len(spt_list)):
		if 'nan' in str(spt_list[i]):
			spt_str.append('Nan')
		else:
			if '4' in str(spt_list[i])[1]:
				spt_str.append('F'+str(spt_list[i])[3:])
			if '5' in str(spt_list[i])[1]:
				spt_str.append('G'+str(spt_list[i])[3:])		
			if '6' in str(spt_list[i])[1]:
				spt_str.append('K'+str(spt_list[i])[3:])
			if '7' in str(spt_list[i])[1]:
				spt_str.append('M'+str(spt_list[i])[3:])
			else:
				spt_str.append('NOTHING')

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
	spt = Column(np.around(np.array(wphottab['spt'][syn_ind[flagtest]]),2),name='SPT')
	sptst = Column(np.array(spt_str)[syn_ind[flagtest]],name='SPT_STR')
	final_tab.add_columns([ra,dec,w,we,j,je,h,he,av,spt,sptst])
	final_tab.write('SerpensSouth_WJHpy_spt.dat',format='ascii',overwrite=True)

	return final_tab



def spt_histogram(final_tab):

    str_spt = ('O','B','A','F','G','K','M','L','T')

    spt = final_tab['SPT_STR']
    spt = [a[0] for a in spt]
    new_spt = []
    #convert string spt to numerical spt 
    for i in range(len(spt)):
        if 'O' in spt[i]:
            new_spt.append(1.0)
        elif 'B' in spt[i]:
            new_spt.append(2.0)
        elif 'A' in spt[i]:
            new_spt.append(3.0)
        elif 'F' in spt[i]:
            new_spt.append(4.0)
        elif 'G' in spt[i]:
            new_spt.append(5.0)
        elif 'K' in spt[i]:
            new_spt.append(6.0)
        elif 'M' in spt[i]:
            new_spt.append(7.0)
        elif 'L' in spt[i]:
            new_spt.append(8.0)
        elif 'T' in spt[i]:
            new_spt.append(9.0)

    bins = np.linspace(1,10,10)
    plt.figure()
    plt.hist(new_spt,bins,facecolor='paleturquoise',edgecolor='k')
    plt.xticks((1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),('O','B','A','F','G','K','M','L','T'))
    plt.xlabel('SpT')
    plt.ylabel('N')
    plt.title('Distribution of Spectral Types')

    plt.show()

    return

# - - - - - - - - 
def main():

	# Read in photometry tables
	phottab = Table.read('SerpensSouth_ALL-PHOT.csv',format='csv')
	wphottab = Table.read('SerpensSouth_W-PHOT_spt.dat',format='ascii')
	w = Table.read('serpenssouth_sexcat_w.fits',format='fits')

	tab = find_zeropoint(phottab,wphottab,w) 
	spt_histogram(tab)
	return tab 

tab = main()
