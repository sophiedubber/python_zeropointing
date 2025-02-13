# - - - - - - - - -
#code to convert 1st part of IDL zeropointing code to python3
# - - - - - - - - 
import sys
import numpy as np
import pandas as pd
import scipy.io as scio
import matplotlib.pyplot as plt

from astropy.io import fits,ascii
from astropy.table import Table,Column
from astropy import units as u
from astropy.coordinates import SkyCoord

# - - - - - - - - -
# Input files: arg[1],arg[2]
# - J: serpsensouth_sexcat_j.fits
# - H: serpenssouth_sexcat_h.fits

# Output files: arg[3]
# - SigmaOri_JH.csv
# - - - - - - - - - 


# - - - - -  F U N C T I O N S - - - - -
#
# FUNCTION to calculate zerp-point std for a specific band
# 
def get_zpstd(a,vec1,uplim,off1,off2,off3):
	# 'Transforms using relations determined from SpeX sythetic photometry, start with J band
	good,offset0 = [],[]
	for i in range(len(a)):
		if a['VECTOR_ASSOC'][i][vec1] < uplim and a['VECTOR_ASSOC'][i][0] > 8. and a['VECTOR_ASSOC'][i][1] >8. and a['FLAGS'][i] == 0 and a['CLASS_STAR'][i] > 0.5: 
			good.append(i)
		jh = (a['VECTOR_ASSOC'][i][0] - a['VECTOR_ASSOC'][i][1])
		offset0.append(off1+off2*(jh) + off3*(jh)**2 + a['VECTOR_ASSOC'][i][vec1] - a['MAG_AUTO'][i])
	# Iterate fit 
	for i in range(5):
		offset0 = np.asarray(offset0)
		zpoff = np.median(offset0[good])
		zperr = np.median(abs(offset0[good]-zpoff))
		good = []
		for k in range(len(a)):
        		if a['VECTOR_ASSOC'][k][vec1] < uplim and a['VECTOR_ASSOC'][k][0] > 8. and a['VECTOR_ASSOC'][k][1] >8. and a['FLAGS'][k] == 0 and a['CLASS_STAR'][k] > 0.5 and abs(offset0[k]-zpoff) < 4.5*zperr: 
                		good.append(k)

	zpstd = np.std(offset0[good])/np.sqrt(float(len(good)-1))
	return zpstd,zpoff

#
# FUNCTION that uses astropy.coordinates to match the two catalgoues
#
def cat_match(c1,c2,max_sep,colra1,coldec1,colra2,coldec2,CONVERT1=True,CONVERT2=True):

	if CONVERT1:
		cat1 = SkyCoord(ra=c1[colra1]*u.degree, dec=c1[coldec1]*u.degree)
	else:
		cat1 = SkyCoord(ra=c1[colra1],dec=c1[coldec1])
	if CONVERT2:
		cat2 = SkyCoord(ra=c2[colra2]*u.degree, dec=c2[coldec2]*u.degree)
	else:
		cat2 = SkyCoord(ra=c2[colra2], dec=c2[coldec2])
	idx, d2d, d3d = cat1.match_to_catalog_sky(cat2)
	sep_constraint = np.where(d2d < max_sep)
	cat1_matches = c1[sep_constraint]
	cat2_matches = c2[idx[sep_constraint]]

	return cat1_matches,cat2_matches,sep_constraint[0],idx[sep_constraint]

# 
# FUNCTION that calculates photometric table columns and saves table
#
def sav_tab(jm,hm,jzpoff,hzpoff,jzpstd,hzpstd):
	
	jmag = [a+jzpoff for a in jm['MAG_AUTO']]
	hmag = [a+hzpoff for a in hm['MAG_AUTO']]
	jmage = [np.sqrt(a**2 + jzpstd**2) for a in jm['MAGERR_AUTO']]	
	hmage = [np.sqrt(a**2 + hzpstd**2) for a in hm['MAGERR_AUTO']]
	rd = [0 for i in range(len(hm))]
	dd = [0 for i in range(len(hm))]

	for i in range(len(jm)):
		rd[i] = np.mean([hm['ALPHA_J2000'][i],jm['ALPHA_J2000'][i]])
		dd[i] = np.mean([hm['DELTA_J2000'][i],jm['DELTA_J2000'][i]])

	js = jm['CLASS_STAR']
	hs = hm['CLASS_STAR']

	sav_dict = {'rd':rd,'dd':dd,'jm':jmag,'jme':jmage,'hm':hmag,'hme':hmage,'js':js,'hs':hs}
	df = pd.DataFrame(sav_dict,columns=['rd','dd','jm','jme','hm','hme','js','hs'])
	df.to_csv(str(sys.argv[3]))

	return

# - - - - - - - - - - - - - - - - - - 

def main():
	# Read in input files
	file1 = str(sys.argv[1])
	file2 = str(sys.argv[2])

	#Check input file extension:
	if '.sexcat' in file1:
		j = Table.read(file1,hdu=2)
		h = Table.read(file2,hdu=2)

	elif '.fits' in file1:
		j = Table.read(file1,format='fits')
		h = Table.read(file2,format='fits')

	jzpstd,jzpoff = get_zpstd(j,0,15.8,-0.0047,-0.0448,0.00964)
	print(jzpstd)
	hzpstd,hzpoff = get_zpstd(h,1,15.1,-0.0056,0.0519,-0.00883)
	print(hzpstd)

	jmed = 2*np.median(j['FWHM_IMAGE'])
	hmed = 2*np.median(h['FWHM_IMAGE'])
	# Get subset of data with good detections
	jgood = np.where((j['FLAGS'] == 0)&(j['FWHM_IMAGE'] > 0.5)&(j['FWHM_IMAGE'] < jmed))
	hgood = np.where((h['FLAGS'] == 0)&(h['FWHM_IMAGE'] > 0.5)&(h['FWHM_IMAGE'] < hmed))

	jcat = j[jgood]
	hcat = h[hgood]
	# Crossmatch J and H catalogues: here the IDL code uses kna_match_cats - could I make use of astropy instead?
	hm,jm,indh,indj = cat_match(hcat,jcat,1.0*u.arcsec,'ALPHA_J2000','DELTA_J2000','ALPHA_J2000','DELTA_J2000',CONVERT1=False,CONVERT2=False)	
	print(len(hcat),' sources in 1')
	print(len(jcat),' sources in 2')
	print(len(hm),' sources matched')


	sav_tab(jm,hm,jzpoff,hzpoff,jzpstd,hzpstd)

	return 

main()

