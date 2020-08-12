# - - - - - - - - 
# Code to convert 2nd part of IDL zeropointing code to python3
# - - - - - - - - 
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import astropy.coordinates as coord

from astroquery.vizier import Vizier
from astropy.table import Table,Column
from make_jh_sav import cat_match
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
	print(catalog)
	tab = catalog[viznam]

	return tab

def matching(matchcat,JH,bands,edit_list):

        if len(matchcat) > 0:
                match,JH_match,ind_match,ind_JH = cat_match(matchcat,JH,mr*u.arcsec,'RAJ2000','DEJ2000','rd','dd',CONVERT=False)
                for i in range(len(edit_list)):
                        bands[edit_list[i][0]][ind_JH] = matchcat[edit_list[i][1]][ind_match]

	return bands

def make_table(JH,tm,wise,usno,apass,sdss,denis,ukidsg,ukidsl,mr):
	
	# Find matches between JH table and 2MASS
	# ind_tm = indices of objects that are matched in 2MASS, ind_JH = indices of objects that are matched in WBand
	tm_match,JH_match,ind_tm,ind_JH = cat_match(tm,JH,mr*u.arcsec,'RAJ2000','DEJ2000','rd','dd',CONVERT=False)
	injh = np.asarray([0 for i in range(len(tm))])
	injh[ind_tm] = 1
	# Unique 2MASS detections - objects in 2MASS query of region that AREN'T in WBand catalogue
	indexu2m = np.where(injh == 0)

	# Extend JH magnitude columns by number of unique 2MASS detections
	b = [np.nan for i in range(len(indexu2m))]
	jm = [*(JH['jm'].tolist()),*b]
	jme = [*(JH['jme'].tolist()),*b]
	hm = [*(JH['hm'].tolist()),*b]
	hme = [*(JH['hme'].tolist()),*b]
	#print(jm)
	a = np.asarray([np.nan for i in range(len(JH))])
	bands = {'j2m':a,'j2me':a,'h2m':a,'h2me':a,'k2m':a,'k2me':a,'ksm':a,'ksme':a,'bm':a,'bme':a,'vm':a,'vme':a,'rcm':a,'rcme':a,'icm':a,'icme':a,'gm':a,'gme':a,'rm':a,'rme':a,'im':a,'ime':a,'zm':a,'zme':a,'zps1m':a,'zps1me':a,'yps1m':a,'yps1me':a,'zum':a,'zume':a,'yum':a,'yume':a,'w1m':a,'w1me':a,'w2m':a,'w2me':a}

	# Store 2MASS magnitude for each object in new columns
	edit_list = [['j2m','Jmag'],['h2m','Hmag'],['k2m','Kmag'],['j2me','e_Jmag'],['h2me','e_Hmag'],['k2me','e_Kmag']]
	for i in range(len(edit_list)):
		bands[edit_list[i][0]][ind_JH] = tm[edit_list[i][1]][ind_tm]
		# Add in 2MASS mags for unique background objects - extend arrays
		bands[edit_list[i][0]] = [*(bands[edit_list[i][0]].tolist()),*(np.asarray(tm['Jmag'])[indexu2m].tolist())]

	# Find matches between WISE and JH
	edit_list = [['w1m','W1mag'],['w1me','e_W1mag'],['w2m','W2mag'],['w2me','e_W2mag']]
	bands = matching(wise,JH,bands,edit_list)

	# Find matches between USNO and JH
	edit_list = [['rcm','R2mag'],['icm','Imag']]
	bands = matching(usno,JH,bands,edit_list)
	bands['rcme'][ind_JH] = 0.3
	bands['icme'][ind_JH] = 0.3

	# Find matches between APASS and JH
	edit_list = [['bm','Bmag'],['bme','e_Bmag'],['vm','Vmag'],['vme','e_Vmag'],['gm','g_mag'],['gme','e_g_mag'],['rm','r_mag'],['rme','e_r_mag'],['im','i_mag'],['ime','e_i_mag']]
	bands = matching(apass,JH,bands,edit_list)

	# Find matches between SDSS and JH
	edit_list = [['gm','gmag'],['gme','e_gmag'],['rm','rmag'],['rme','e_rmag'],['im','imag'],['ime','e_imag'],['zm','zmag'],['zme','e_zmag']]
	bands = matching(sdss,JH,bands,edit_list)
	
	# Find matches between DENIS and JH
	edit_list = [['icm','imag'],['icme','e_imag']]
	bands = matching(denis,JH,bands,edit_list)
	finite_check = np.where(np.isfinite(bands['k2m'][ind_JH]))
	return 


def main():
	# Read in JH catalogue
	JH = Table.read('SerpensSouth_JH.csv',format='csv')
	# Matching radius (arcsecs)
	mr = 2.0
	# Option to just look at core of region
	core = np.where((JH['rd'] > 0.)&(JH['rd']<360.))
	cen = [np.mean([max(JH['rd'][core]),min(JH['rd'][core])]),np.mean([max(JH['dd'][core]),min(JH['dd'][core])])]
	dis = [((max(JH['rd'][core])-min(JH['rd'][core]))/0.98)*np.cos((cen[1])*(np.pi/180.))*60.,(max(JH['dd'][core])-min(JH['dd'][core]))/0.98*60.]
	print('FIND 2MASS')
	tm = query(cen,dis,'2MASS-PSC','II/246/out')
	wise = query(cen,dis,'allwise','II/328/allwise')
	usno = query(cen,dis,'USNO-B1','I/284/out')
	apass = query(cen,dis,'apass9','II/336/apass9')
	#sdss = query(cen,dis,'V/139/sdss9','V/139/sdss9')
	denis = query(cen,dis,'B/denis','B/denis/denis')
	ukidsg = query(cen,dis,'II/319/gcs9','II/319/gcs9')
	ukidsl = query(cen,dis,'II/319/las9','II/319/las9')
	make_table(JH,tm,mr)

	return

main()
