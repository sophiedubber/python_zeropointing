# - - - - - - - - -
#code to convert 1st part of IDL zeropointing code to python3
# - - - - - - - - 
import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt

from astropy.io import fits,asciis
from astropy.table import Table,Column

# - - - - - - - - -
# Input files:
# - J: serpsensouth_sexcat_j.fits
# - H: serpenssouth_sexcat_h.fits

# Output files: 
# - serpenssouth_JH.sav
# - - - - - - - - - 

# Read in input files
j = Table.read('serpenssouth_sexcat_j.fits',format='fits')
h = Table.read('serpenssouth_sexcat_h.fits',format='fits')

# 'Transforms using relations determined from SpeX sythetic photometry, start with J band
j_good,joffset0 = [],[]
for i in range(len(j)):
	if j['VECTOR_ASSOC'][i][0] < 15.8 and j['VECTOR_ASSOC'][i][0] > 8. and j['VECTOR_ASSOC'][i][1] >8. and j['FLAGS'][i] == 0 and j['CLASS_STAR'][i] > 0.5: 
		j_good.append(i)
	jh = (j['VECTOR_ASSOC'][i][0] - j['VECTOR_ASSOC'][i][1])
	joffset0.append(-0.0047-0.0448*(jh) + 0.00964*(jh)**2 + j['VECTOR_ASSOC'][i][0] - j['MAG_AUTO'][i])

# Iterate fit 
for i in range(5):
	joffset0 = np.asarray(joffset0)
	jzpoff = np.median(joffset0[j_good])
	jzperr = np.median(abs(joffset0[j_good]-jzpoff))
	print(jzpoff,jzperr,len(j_good))
	j_good = []
	for k in range(len(j)):
        	if j['VECTOR_ASSOC'][k][0] < 15.8 and j['VECTOR_ASSOC'][k][0] > 8. and j['VECTOR_ASSOC'][k][1] >8. and j['FLAGS'][k] == 0 and j['CLASS_STAR'][k] > 0.5 and abs(joffset0[k]-jzpoff) < 4.5*jzperr: 
                	j_good.append(k)

jzpstd = np.std(joffset0[j_good])/np.sqrt(float(len(j_good)-1))
print(jzpstd)
