# - - - - - - - - 
# Code to convert 2nd part of IDL zeropointing code to python3
# - - - - - - - - 
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import astropy.coordinates as coord

from astroquery.vizier import Vizier
from astropy.table import Table,Column
from make_jh_sav.py import cat_match
Vizier.ROW_LIMIT = 99999
# - - - - - - - - 
# Input Files:
# - SerpensSouth_JH.csv

# Output Files:
# - SerpensSouth_ALL-PHOT.csv
# - SerpensSouth_W-PHOT.csv
# - - - - - - - -

# - - - - - - - - F U N C T I O N S - - - - - - - -
#
# FUNCTION to query each catalogue and create full table of literature photometry
#
def query(cen,dis,cat,viznam):
	# First, query 2mass via vizier
	catalog = Vizier.query_region(coord.SkyCoord(ra=cen[0],dec=cen[1],unit=(u.deg,u.deg),frame='icrs'),width=dis[0]*u.arcmin,height=dis[1]*u.arcmin,catalog=cat)
	tab = catalog[viznam]

	return tab

def make_table(JH,tm,mr):
	# Find matches between JH table and 2MASS
	tm_match,JH_match,ind_tm,ind_JH = cat_match(tm,JH,mr,'raj2000','dej2000','rd','dd')
	injh = [int(0) for i in range(len(tm_match))]
	injh[ind_tm] = 1
	# Unique 2MASS detections
	indexu2m = np.where(injh == 0)

	b = [np.nan for i in range(len(indexu2m))]
	jm = [*(JH['jm'].tolist()),*b]
	jme = [*(JH['jme'].tolist()),*b]
	hm = [*(JH['hm'].tolist()),*b]
	hme = [*(JH['hme'].tolist()),*b]

	a = [np.nan for i in range(len(JH))]
	bands = {'j2m':a,'j2me':a,'h2m':a,'h2me':a,'k2m':a,'k2me':a,'ksm':a,'ksme':a,'bm':a,'bme':a,'vm':a,'vme':a,'rcm':a,'rcme':a,'icm':a,'icme':a,'gm':a,'gme':a,'rm':a,'rme':a,'im':a,'ime':a,'zm':a,'zme':a,'zps1m':a,'zps1me':a,'yps1m':a,'yps1me':a,'zum':a,'zume':a,'yum':a,'yume':a,'w1m':a,'w1me':a,'w2m':a,'w2me':a}

	edit_list = [['j2m','jmag'],['h2m','hmag'],['k2m','kmag'],['j2me','e_jmag'],['h2me','e_hmag'],['k2me','e_kmag']]
	for i in range(len(edit_list)):
		bands[edit_list[i][0]][indexjh] = tm[edit_list[i][1]][indexcat2m]

	# Up to line 48 in IDL


def main():
	# Read in JH catalogue
	JH = Table.read('SerpensSouth_JH.csv',format='csv')
	# Matching radius (arcsecs)
	mr = 2.0
	# Option to just look at core of region
	core = np.where((JH['rd'] > 0.)&(JH['rd']<360.))
	cen = [np.mean([max(JH['rd'][core]),min(JH['rd'][core])]),np.mean([max(JH['dd'][core]),min(JH['dd'][core])])]
	dis = [((max(JH['rd'][core])-min(JH['rd'][core]))/0.98)*np.cos((cen[1])*(np.pi/180.))*60.,(max(JH['dd'][core])-min(JH['dd'][core]))/0.98*60.]

	tm = query(cen,dis,'2MASS-PSC','II/246/out')
	make_table(JH,tm)


	return tmass

t=main()
