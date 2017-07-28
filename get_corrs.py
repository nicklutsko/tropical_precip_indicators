"""

Script for calculating correlations in the indicators from CMIP, AMIP and AQUA simulations

Nick Lutsko and Max Popp 2017

"""
from scipy.io import netcdf
import numpy as np
import scipy.stats as ss

print "Correlations between indicators"

#Different sets of simulations
opts = [ "CMIP Historical 1979 - 2005", "CMIP Historical 1861 - 2005", "AMIP 1979 - 2005", "AMIP 1979 - 2008", "aqua 1979 - 1983 " ]

#File names
files = [ 'Precipitation_indicators_SAH_historical_r1i1p1_197901-200512_ZM.nc', 'Precipitation_indicators_SAH_historical_r1i1p1_186101-200512_ZM.nc', 'Precipitation_indicators_SAH_amip_r1i1p1_197901-200512_ZM.nc', 
'Precipitation_indicators_SAH_amip_r1i1p1_197901-200812_ZM.nc',
'Precipitation_indicators_SAH_aqua_r1i1p1_197901-198312_ZM.nc' ]

#Indicators
let = [ 'E', 'D', 'A', 'H', 'R', 'I', 'S', 'L', 'W', 'Ep', 'Ap', 'Mp' ]
l = len( let )

for z in range( 5 ):

	print opts[z]

	#Load indicator file
	f1 = netcdf.netcdf_file( files[z], 'r' )

	if z < 2: #CMIP
		inds = np.zeros( ( l, 26 ) )
	elif z > 1 and z < 4: #AMIP
		inds = np.zeros( ( l, 24 ) )
	else: #AQUA
		inds = np.zeros( ( l, 8 ) )

	#Load indicators
	for i in range( l ):
		inds[i] = f1.variables[ let[i] ][:]

	#Print correlations, alerting if correlation is statistically significant
	for i in range( l ):
		for j in range( l ):
			if ss.linregress( inds[i], inds[j] )[3] < 0.025:
				print let[i], let[j], "%.2f" % ss.linregress( inds[i], inds[j] )[2] ** 2, "sig"
			else:
				print let[i], let[j], "%.2f" % ss.linregress( inds[i], inds[j] )[2] ** 2


