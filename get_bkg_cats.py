# - - - - - - - - 
# Code to convert 2nd part of IDL zeropointing code to python3
# - - - - - - - - 
import numpy as np
import pandas as pd
import astropy.units as u
import scipy.io as scio
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
	if len(catalog) > 0:
		tab = catalog[viznam]
		print(str(len(tab))+' sources matched')
	else: 	
		tab = Table()
		print('No sources matched')
	return tab

#
# FUNCTION to add in bands from catalogues to main photometry table
#
def matching(matchcat,JH,mr,bands,edit_list):

	if len(matchcat) > 0:
		match,JH_match,ind_match,ind_JH = cat_match(matchcat,JH,mr*u.arcsec,'RAJ2000','DEJ2000','rd','dd',CONVERT=False)
		for i in range(len(edit_list)):
			bands[edit_list[i][0]][ind_JH] = matchcat[edit_list[i][1]][ind_match]

	return bands

#
# FUNCTION to read in and store data from large .sav table containing IRTF observed and synthetic photometry 
def IRTF_tab_read():

	a = scio.readsav('../IRTF_SpeX_phot.sav',python_dict=True)
	b = a['spexphot']
	c = a['obsphot']
	irtf_spexphot = Table(b)
	irtf_obsphot = Table(c)

	return irtf_spexphot,irtf_obsphot
	
#
# FUNCTION to crossmatch with 2 mass and then add all photometric bands into one big table
#
def make_table(JH,tm,wise,usno,apass,sdss,denis,ukidsg,ukidsl,panstars,mr):
	
	# Find matches between JH table and 2MASS
	# ind_tm = indices of objects that are matched in 2MASS, ind_JH = indices of objects that are matched in WBand
	tm_match,JH_match,ind_tm,ind_JH = cat_match(tm,JH,mr*u.arcsec,'RAJ2000','DEJ2000','rd','dd',CONVERT=False)
	injh = np.asarray([0 for i in range(len(tm))])
	injh[ind_tm] = 1
	# Unique 2MASS detections - objects in 2MASS query of region that AREN'T in WBand catalogue
	indexu2m = np.where(injh == 0)
	# Extend JH magnitude columns by number of unique 2MASS detections
	b = [np.nan for i in range(len(indexu2m[0]))]
	jm = [*(JH['jm'].tolist()),*b]
	jme = [*(JH['jme'].tolist()),*b]
	hm = [*(JH['hm'].tolist()),*b]
	hme = [*(JH['hme'].tolist()),*b]

	a = np.asarray([np.nan for i in range(len(jm))])
	b = np.asarray([np.nan for i in range(len(JH))])
	bands = {'rd':JH['rd'],'dd':JH['dd'],'jm':jm,'jme':jme,'hm':hm,'hme':hme,'j2m':b,'j2me':b,'h2m':b,'h2me':b,'k2m':b,'k2me':b,'ksm':a,'ksme':a,'bm':a,'bme':a,'vm':a,'vme':a,'rcm':a,'rcme':a,'icm':a,'icme':a,'gm':a,'gme':a,'rm':a,'rme':a,'im':a,'ime':a,'zm':a,'zme':a,'zps1m':a,'zps1me':a,'yps1m':a,'yps1me':a,'zum':a,'zume':a,'yum':a,'yume':a,'w1m':a,'w1me':a,'w2m':a,'w2me':a}

	# Extend coordinate arrays
	bands['rd'] = [*(bands['rd'].tolist()),*(np.asarray(tm['RAJ2000'])[indexu2m].tolist())]
	bands['dd'] = [*(bands['dd'].tolist()),*(np.asarray(tm['DEJ2000'])[indexu2m].tolist())]

	# Store 2MASS magnitude for each object in new columns
	edit_list = [['j2m','Jmag'],['h2m','Hmag'],['k2m','Kmag'],['j2me','e_Jmag'],['h2me','e_Hmag'],['k2me','e_Kmag']]
	for i in range(len(edit_list)):
		bands[edit_list[i][0]][ind_JH] = tm[edit_list[i][1]][ind_tm]
		# Add in 2MASS mags for unique background objects - extend arrays
		bands[edit_list[i][0]] = np.asarray([*(bands[edit_list[i][0]].tolist()),*(np.asarray(tm['Jmag'])[indexu2m].tolist())])
	# Make new coordinate table for matching
	JH = Table()
	racol = Column(bands['rd'],name='rd')
	decol = Column(bands['dd'],name='dd')
	JH.add_columns([racol,decol])

	# Find matches between WISE and JH
	edit_list = [['w1m','W1mag'],['w1me','e_W1mag'],['w2m','W2mag'],['w2me','e_W2mag']]
	bands = matching(wise,JH,mr,bands,edit_list)

	# Find matches between USNO and JH
	edit_list = [['rcm','R2mag'],['icm','Imag']]
	bands = matching(usno,JH,mr,bands,edit_list)
	bands['rcme'][ind_JH] = 0.3
	bands['icme'][ind_JH] = 0.3

	# Find matches between APASS and JH
	edit_list = [['bm','Bmag'],['bme','e_Bmag'],['vm','Vmag'],['vme','e_Vmag'],['gm','g_mag'],['gme','e_g_mag'],['rm','r_mag'],['rme','e_r_mag'],['im','i_mag'],['ime','e_i_mag']]
	bands = matching(apass,JH,mr,bands,edit_list)

	# Find matches between SDSS and JH
	edit_list = [['gm','gmag'],['gme','e_gmag'],['rm','rmag'],['rme','e_rmag'],['im','imag'],['ime','e_imag'],['zm','zmag'],['zme','e_zmag']]
	bands = matching(sdss,JH,mr,bands,edit_list)
	
	# Find matches between DENIS and JH
	edit_list = [['icm','Imag'],['icme','e_Imag']]
	if len(denis) > 0:
		denis_match,JH_match,ind_denis,ind_JH = cat_match(denis,JH,mr*u.arcsec,'RAJ2000','DEJ2000','rd','dd',CONVERT=False)
		for i in range(len(edit_list)):
			bands[edit_list[i][0]][ind_JH] = denis[edit_list[i][1]][ind_denis]
	finite_check = np.where(np.isfinite(np.asarray(bands['k2m'])[ind_JH]))
	if len(finite_check) >= 1:
		bands['k2m'][ind_JH[finite_check]] = denis['Kmag'][ind_denis[finite_check]]
		bands['k2me'][ind_JH[finite_check]] = denis['e_Kmag'][ind_denis[finite_check]]

	# Find matches between UKIDSS_GCS and JH
	edit_list = [['zum','Zmag'],['zume','e_Zmag'],['yum','Ymag'],['yume','e_Ymag']]
	bands = matching(ukidsg,JH,mr,bands,edit_list)

	# Find matches between UKIDSS_LAS and JH
	edit_list = [['yum','Ymag'],['yume','e_Ymag']]
	bands = matching(ukidsl,JH,mr,bands,edit_list)

	# Find matches between PANSTARRS-1 and JH
	edit_list = [['gm','gmag'],['gme','e_gmag'],['rm','rmag'],['rme','e_rmag'],['im','imag'],['ime','e_imag'],['zps1m','zmag'],['zps1me','e_zmag'],['yps1m','ymag'],['yps1me','e_ymag']]
	bands = matching(panstars,JH,mr,bands,edit_list)

	# Save table with all photometry
	df = pd.DataFrame(bands,columns=['rd','dd','jm','jme','hm','hme','j2m','j2me','h2m','h2me','k2m','k2me','ksm','ksme','bm','bme','vm','vme','rcm','rcme','icm','icme','gm','gme','rm','rme','im','ime','zm','zme','zps1m','zps1me','yps1m','yps1me','zum','zume','yum','yume','w1m','w1me','w2m','w2me'])
	#df.to_csv('SerpensSouth_ALL-PHOT.csv')

	return bands

u
# FUNCTION to calculate W photometry using full photometric table
#
def get_w(bands):

	spexphot,obsphot = IRTF_tab_read()
	# Create full arrays of photometry and errors
	wobsphot = [[bands['jm']],[bands['hm']],[bands['ksm']],[bands['j2m']],[bands['h2m']],[bands['k2m']],[bands['w1m']],[bands['w2m']],[bands['bm']],[bands['vm']],[bands['rcm']],[bands['icm']],[bands['gm']],[bands['rm']],[bands['im']],[bands['zm']],[bands['zps1m']],[bands['yps1m']],[bands['zum']],[bands['yum']]]
	wobsphot = [[bands['jme']],[bands['hme']],[bands['ksme']],[bands['j2me']],[bands['h2me']],[bands['k2me']],[bands['w1me']],[bands['w2me']],[bands['bme']],[bands['vme']],[bands['rcme']],[bands['icme']],[bands['gme']],[bands['rme']],[bands['ime']],[bands['zme']],[bands['zps1me']],[bands['yps1me']],[bands['zume']],[bands['yume']]]
	stdphot = [[spexphot['W']],[spexphot['WC']],[spexphot['JMKO']],[spexphot['HMKO']],[spexphot['KSMKO']],[spexphot['J2M']],[spexphot['H2M']],[spexphot['K2M']],[obsphot['W1']],[obsphot['W2']],[obsphot['BJ']],[obsphot['VJ']],[obsphot['RC']],[obsphot['IC']],[obsphot['G_SD']],[obsphot['R_SD']],[obsphot['I_SD']],[spexphot['ZSD']],[spexphot['ZPS1']],[spexphot['YPS1']],[spexphot['ZUK']],[spexphot['YUK']]] 
	nfilt = len(stdphot)
	# Only use objects with enough photometric bands and av info
#	goodstd = np.where(

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

	# Use astroquery to find WBand objects in literature photometric catalogues
	tm = query(cen,dis,'2MASS-PSC','II/246/out')
	wise = query(cen,dis,'allwise','II/328/allwise')
	usno = query(cen,dis,'USNO-B1','I/284/out')
	apass = query(cen,dis,'apass9','II/336/apass9')
	sdss = query(cen,dis,'SDSS','V/139/sdss9')
	denis = query(cen,dis,'B/denis','B/denis/denis')
	ukidsg = query(cen,dis,'II/319/gcs9','II/319/gcs9')
	ukidsl = query(cen,dis,'II/319/las9','II/319/las9')
	panstars = query(cen,dis,'PS1','II/349/ps1')

	# Make full photometry table
	bands = make_table(JH,tm,wise,usno,apass,sdss,denis,ukidsg,ukidsl,panstars,mr)

	# Get W-Band magnitudes
	get_w(bands)

	return bands

bands = main()
