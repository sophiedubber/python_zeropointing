# - - - - - - - - 
# Code to convert 2nd part of IDL zeropointing code to python3
# - - - - - - - -
import sys
import time 
import numpy as np
import pandas as pd
import astropy.units as u
import scipy.io as scio
import scipy.stats as ss
import matplotlib.pyplot as plt
import astropy.coordinates as coord

from multiprocessing import Pool
from astroquery.vizier import Vizier
from astropy.table import Table,Column,join,vstack
from make_jh_sav import cat_match
from itertools import zip_longest
Vizier.ROW_LIMIT = 99999
# - - - - - - - - 
# Input Files e.g: arg[1]
# - SerpensSouth_JH.csv

# Output Files e.g: arg[2], arg[3]
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
def matching(matchcat,JH,mr,bands,edit_list,ra='RAJ2000',dec='DEJ2000'):

	if len(matchcat) > 0:
		match,JH_match,ind_match,ind_JH = cat_match(matchcat,JH,mr*u.arcsec,ra,dec,'rd','dd',CONVERT1=False)
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
	tm_match,JH_match,ind_tm,ind_JH = cat_match(tm,JH,mr*u.arcsec,'RAJ2000','DEJ2000','rd','dd',CONVERT1=False)
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

	bands = {'rd':JH['rd'],'dd':JH['dd'],'jm':jm,'jme':jme,'hm':hm,'hme':hme,'j2m':np.asarray([np.nan for i in range(len(JH))]),'j2me':np.asarray([np.nan for i in range(len(JH))]),'h2m':np.asarray([np.nan for i in range(len(JH))]),'h2me':np.asarray([np.nan for i in range(len(JH))]),'k2m':np.asarray([np.nan for i in range(len(JH))]),'k2me':np.asarray([np.nan for i in range(len(JH))]),'ksm':np.asarray([np.nan for i in range(len(jm))]),'ksme':np.asarray([np.nan for i in range(len(jm))]),'bm':np.asarray([np.nan for i in range(len(jm))]),'bme':np.asarray([np.nan for i in range(len(jm))]),'vm':np.asarray([np.nan for i in range(len(jm))]),'vme':np.asarray([np.nan for i in range(len(jm))]),'rcm':np.asarray([np.nan for i in range(len(jm))]),'rcme':np.asarray([np.nan for i in range(len(jm))]),'icm':np.asarray([np.nan for i in range(len(jm))]),'icme':np.asarray([np.nan for i in range(len(jm))]),'gm':np.asarray([np.nan for i in range(len(jm))]),'gme':np.asarray([np.nan for i in range(len(jm))]),'rm':np.asarray([np.nan for i in range(len(jm))]),'rme':np.asarray([np.nan for i in range(len(jm))]),'im':np.asarray([np.nan for i in range(len(jm))]),'ime':np.asarray([np.nan for i in range(len(jm))]),'zm':np.asarray([np.nan for i in range(len(jm))]),'zme':np.asarray([np.nan for i in range(len(jm))]),'zps1m':np.asarray([np.nan for i in range(len(jm))]),'zps1me':np.asarray([np.nan for i in range(len(jm))]),'yps1m':np.asarray([np.nan for i in range(len(jm))]),'yps1me':np.asarray([np.nan for i in range(len(jm))]),'zum':np.asarray([np.nan for i in range(len(jm))]),'zume':np.asarray([np.nan for i in range(len(jm))]),'yum':np.asarray([np.nan for i in range(len(jm))]),'yume':np.asarray([np.nan for i in range(len(jm))]),'w1m':np.asarray([np.nan for i in range(len(jm))]),'w1me':np.asarray([np.nan for i in range(len(jm))]),'w2m':np.asarray([np.nan for i in range(len(jm))]),'w2me':np.asarray([np.nan for i in range(len(jm))])}

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
	bands = matching(sdss,JH,mr,bands,edit_list,ra='RA_ICRS',dec='DE_ICRS')
	
	# Find matches between DENIS and JH
	edit_list = [['icm','Imag'],['icme','e_Imag']]
	if len(denis) > 0:
		denis_match,JH_match,ind_denis,ind_JH = cat_match(denis,JH,mr*u.arcsec,'RAJ2000','DEJ2000','rd','dd',CONVERT1=False)
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
	df.to_csv(str(sys.argv[2]))

	return bands

#
# FUNCTION to calculate W photometry using full photometric table
#
def get_w(input_data):

	wobsphot = input_data[0]
	wobsphote = input_data[1]

	spexphot,obsphot = IRTF_tab_read()
	stdphot = np.asarray([spexphot['W'],spexphot['WC'],spexphot['JMKO'],spexphot['HMKO'],spexphot['KSMKO'],spexphot['J2M'],spexphot['H2M'],spexphot['K2M'],obsphot['W1'],obsphot['W2'],obsphot['BJ'],obsphot['VJ'],obsphot['RC'],obsphot['IC'],obsphot['G_SD'],obsphot['R_SD'],obsphot['I_SD'],spexphot['ZSD'],spexphot['ZPS1'],spexphot['YPS1'],spexphot['ZUK'],spexphot['YUK'],spexphot['INUM']])
	stdphot_rows = list([x for x in y if x is not None] for y in zip_longest(*stdphot))
	goodstd = np.where(np.logical_and(np.asarray([sum(np.isfinite(a)) for a in stdphot_rows]) > 5,np.logical_or(np.asarray(spexphot['AV']) > 0,np.asarray(spexphot['RV']) == 3.1)))
	gspexphot = spexphot[goodstd[0]]	
	gobsphot = obsphot[goodstd[0]]
	gstdphot_rows = np.asarray(stdphot_rows)[goodstd[0]]

	gstdphot = np.asarray([gspexphot['W'],gspexphot['WC'],gspexphot['JMKO'],gspexphot['HMKO'],gspexphot['KSMKO'],gspexphot['J2M'],gspexphot['H2M'],gspexphot['K2M'],gobsphot['W1'],gobsphot['W2'],gobsphot['BJ'],gobsphot['VJ'],gobsphot['RC'],gobsphot['IC'],gobsphot['G_SD'],gobsphot['R_SD'],gobsphot['I_SD'],gspexphot['ZSD'],gspexphot['ZPS1'],gspexphot['YPS1'],gspexphot['ZUK'],gspexphot['YUK'],gspexphot['INUM']])

	njh = len(wobsphot[0])
	nstd = len(gstdphot[0])
	nfilt = len(gstdphot)
	results_tab=Table()
	west=[np.nan for i in range(njh)]
	weste=[np.nan for i in range(njh)]
	wcest=[np.nan for i in range(njh)]
	wceste=[np.nan for i in range(njh)]
	avest=[np.nan for i in range(njh)]
	rd = [np.nan for i in range(njh)]
	dd = [np.nan for i in range(njh)]
	spt = [np.nan for i in range(njh)]
	for i in range(njh):
		if sum(np.isfinite(wobsphot[6:,i])) >= 1:
			# Normalise to J H K
			tab1 = np.asarray([[a for j in range(nstd)] for a in wobsphot[2:7,i]])
			norm = np.nanmean(tab1-gstdphot[2:7],axis=0)
			tab2 = np.asarray([[a for j in range(nfilt)] for a in norm])
			tab2_rows = list([x for x in y if x is not None] for y in zip_longest(*tab2))
			normphot  = gstdphot + np.asarray(tab2_rows)

			# Simple MCMC to get expected W
			temperrphot = np.choose(wobsphote[2:,i] < 0.02, (wobsphote[2:,i], 0.02))
			nsims = 20
			temp_w,temp_wc,temp_av,temp_spt = [],[],[],[]
			for k in range(nsims):
				tempophot = wobsphot[2:,i] + np.random.normal(size=nfilt-3)*temperrphot
				tab1 = np.asarray([[a for j in range(nstd)] for a in tempophot])
				tab2 = np.asarray([[a for j in range(nstd)] for a in temperrphot])
				offset = (normphot[2:nfilt-1,:] - tab1)/tab2
				bestsub = np.argmin(np.nansum(offset**2,axis=0)/np.nansum(np.isfinite(offset),axis=0),axis=0)
				
				temp_w.append(normphot[0,bestsub])
				temp_wc.append(normphot[1,bestsub])
				temp_av.append(gspexphot['AV'][bestsub])
				temp_spt.append(gspexphot['INUM'][bestsub])
				
			west[i] = np.nanmedian(temp_w)
			weste[i] = np.nanstd(temp_w)	
			wcest[i] = np.nanmedian(temp_wc)
			wceste[i] = np.nanstd(temp_wc)
			avest[i] = np.nanmean(temp_av)
			spt[i] = np.float(ss.mode(temp_spt)[0][0])
		rd[i] = wobsphot[0,i]
		dd[i] = wobsphot[1,i]

	wcol = Column(np.around(west,6),name='west')
	wecol = Column(np.around(weste,6),name='weste')
	wccol = Column(np.around(wcest,6),name='wcest')
	wcecol = Column(np.around(wceste,6),name='wceste')
	acol = Column(np.around(avest,6),name='avest')
	rcol = Column(np.around(rd,4),name='rd')
	dcol = Column(np.around(dd,4),name='dd')
	scol = Column(np.around(spt,2),name='spt')

	results_tab.add_columns([rcol,dcol,wcol,wecol,wccol,wcecol,acol,scol])

	return results_tab


# Try to paralellise function that calculates W - previously have used Pool, which requires one input split into equal chunks  - can i supply multiple inputs all split into the same sized chunks?

def mp_get_w(bands):

        print('Reading in IRTF table...')
        spexphot,obsphot = IRTF_tab_read()

	# Create full arrays of photometry and errors
        wobsphot = np.asarray([bands['rd'],bands['dd'],bands['jm'],bands['hm'],bands['ksm'],bands['j2m'],bands['h2m'],bands['k2m'],bands['w1m'],bands['w2m'],bands['bm'],bands['vm'],bands['rcm'],bands['icm'],bands['gm'],bands['rm'],bands['im'],bands['zm'],bands['zps1m'],bands['yps1m'],bands['zum'],bands['yum']])
        wobsphote = np.asarray([bands['rd'],bands['dd'],bands['jme'],bands['hme'],bands['ksme'],bands['j2me'],bands['h2me'],bands['k2me'],bands['w1me'],bands['w2me'],bands['bme'],bands['vme'],bands['rcme'],bands['icme'],bands['gme'],bands['rme'],bands['ime'],bands['zme'],bands['zps1me'],bands['yps1me'],bands['zume'],bands['yume']])
       
        # Setup multiprocessing function: split tables into chunks
        size1 = int(len(wobsphot[0])/4)
        size2 = int(len(wobsphote[0])/4)  

        print('Calculating W...')
        start = time.time()
        data_chunks = [[wobsphot[:,0:size1],wobsphote[:,0:size2]],[wobsphot[:,size1:2*size1],wobsphote[:,size2:2*size2]],[wobsphot[:,2*size1:3*size1],wobsphote[:,2*size2:3*size2]],[wobsphot[:,3*size1:],wobsphote[:,3*size2:]]]
        pool = Pool(processes=4) 
        result = pool.map(get_w,data_chunks)

        # Join results table back together
        results = vstack([result[0],result[1],result[2],result[3]])
        results.write(str(sys.argv[3]),format='ascii',overwrite=True)

        end = time.time()
        print('Took',str(end-start))
        return

# - - - - - - - - 

def main():
	# Read in JH catalogue
	filen = str(sys.argv[1])
	JH = Table.read(filen,format='csv')
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

	# Get W-Band magnitudes and write to file
	sav_dict = mp_get_w(bands)

	return

main()
