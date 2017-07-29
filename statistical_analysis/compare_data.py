"""

Script for comparing correlations of indicators between CMIP, AMIP and AQUA simulations

Nick Lutsko and Max Popp 2017

"""
from scipy.io import netcdf
import numpy as np
import scipy.stats as ss

print "Compare datasets"

#Load files with indicators
files = [ 'Precipitation_indicators_SAH_historical_r1i1p1_197901-200512_ZM.nc', 'Precipitation_indicators_SAH_amip_r1i1p1_197901-200512_ZM.nc', 
'Precipitation_indicators_SAH_aqua_r1i1p1_197901-198312_ZM.nc' ]

f1 = netcdf.netcdf_file( files[0], 'r' )
f2 = netcdf.netcdf_file( files[1], 'r' )
f3 = netcdf.netcdf_file( files[2], 'r' )

#Indicators
let = [ 'E', 'D', 'A', 'H', 'R', 'I', 'S', 'L', 'W', 'Ep', 'Ap', 'Mp' ]
l = len( let )

print "CMIP and AMIP"

#Models with both AMIP and CMIP simulations
o = [0, 1, 2, 3, 4, 5, 7, 9, 10, 13, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25] #CMIP
o2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22, 23] #AMIP

inds = np.zeros( ( ( 2, l, 20 ) ) )

#Get indicators for models with both AMIP and CMIP simulations
a = 0
b = 0
for i in range( 26 ):

	if a < 20:
		#CMIP
		if i == o[a] and a < 20:
			for j in range( l ):
				inds[0, j, a] = f1.variables[ let[j] ][i]
			a += 1
	if b < 20:
		#AMIP
		if i == o2[b] and b < 20:
			for j in range( l ):
				inds[1, j, b] = f2.variables[ let[j] ][i]
			b += 1

#Print correlations
for i in range( l ):
	print let[i], np.corrcoef( inds[0, i], inds[1, i] )[0, 1] ** 2, ss.linregress( inds[0, i], inds[1, i] )[3]

print "CMIP and AQUA"

#Models with AMIP, CMIP and AQUA simulations
o = [ 5, 7, 9, 17, 20, 22, 23, 24]
o2 = [ 5, 6, 7, 14, 17, 18, 19, 22]
o3 = [0, 1, 2, 3, 4, 5, 6, 7 ]

inds = np.zeros( ( ( 3, l, 8 ) ) )

#Get indicators for models with AMIP, CMIP and AQUA simulations
a = 0
b = 0
c = 0
for i in range( 26 ):

	if a < 8:
		if i == o[a]:
			#CMIP
			for j in range( l ):
				inds[0, j, a] = f1.variables[ let[j] ][i]
			a += 1
	if b < 8:
		if i == o2[b]:
			#AMIP
			for j in range( l ):
				inds[1, j, b] = f2.variables[ let[j] ][i]
			b += 1
	if c < 8:
		if i == o3[c]:
			#AQUA
			for j in range( l ):
				inds[2, j, c] = f3.variables[ let[j] ][i]
			c += 1

#Print correlations
for i in range( l ):
	print let[i], np.corrcoef( inds[0, i], inds[2, i] )[0, 1] ** 2, ss.linregress( inds[0, i], inds[2, i] )[3]

print "AMIP and AQUA"
for i in range( l ):
	print let[i], np.corrcoef( inds[2, i], inds[1, i] )[0, 1] ** 2, ss.linregress( inds[2, i], inds[1, i] )[3]







