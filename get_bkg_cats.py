# - - - - - - - - 
# Code to convert 2nd part of IDL zeropointing code to python3
# - - - - - - - - 
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import astropy.coordinates as coord

from astroquery.vizier import Vizier
from astropy.table import Table,Column
# - - - - - - - - 
# Input Files:
# --> SerpensSouth_JH.csv
#
# Output Files:
# --> SerpensSouth_ALL-PHOT.csv
# --> SerpensSouth_W-PHOT.csv
# - - - - - - - -

# - - - - - - - - F U N C T I O N S - - - - - - - -
#
# FUNCTION to query each catalogue and create full table of literature photometry
#
def query(JH,cen,dis):
	# First, query 2mass via vizier
	tmass = Vizier.query_region(coord.SkyCoord(ra=cen[0],dec=cen[1],unit=(u.deg,u.deg),frame='icrs'),width=dis[0]*u.arcmin,height=dis[1]*u.arcmin,catalog='2MASS-PSC')
	return tmass

def main():
	# read in JH catalogue
	JH = Table.read('SerpensSouth_JH.csv',format='csv')
	# matching radius (arcsecs)
	mr = 2.0
	# option to just look at core of region
	core = np.where((JH['rd'] > 0.)&(JH['rd']<360.))
	cen = [np.mean([max(JH['rd'][core]),min(JH['rd'][core])]),np.mean([max(JH['dd'][core]),min(JH['dd'][core])])]
	print(cen)
	dis = [((max(JH['rd'][core])-min(JH['rd'][core]))/0.98)*np.cos((cen[1])*(np.pi/180.))*60.,(max(JH['dd'][core])-min(JH['dd'][core]))/0.98*60.]
	print(dis)
	tmass = query(JH,cen,dis)
	print(tmass)
	return

main()
