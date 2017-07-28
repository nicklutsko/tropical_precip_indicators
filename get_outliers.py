"""

Script for checking for outliers in the indicators from CMIP, AMIP and AQUA simulations

Nick Lutsko and Max Popp 2017

"""
from scipy.io import netcdf
import numpy as np
import scipy.stats as ss

print "Check for outliers"

#Different sets of simulations
opts = [ "CMIP Historical 1979 - 2005", "CMIP Historical 1861 - 2005", "AMIP 1979 - 2005", "AMIP 1979 - 2008", "aqua 1979 - 1983 " ]

#File names
files = [ 'Precipitation_indicators_SAH_historical_r1i1p1_197901-200512_ZM.nc', 'Precipitation_indicators_SAH_historical_r1i1p1_186101-200512_ZM.nc', 'Precipitation_indicators_SAH_amip_r1i1p1_197901-200512_ZM.nc', 
'Precipitation_indicators_SAH_amip_r1i1p1_197901-200812_ZM.nc',
'Precipitation_indicators_SAH_aqua_r1i1p1_197901-198312_ZM.nc' ]

#Indicators
let = [ 'E', 'D', 'A', 'H', 'R', 'I', 'S', 'L', 'W', 'Ep', 'Ap', 'Mp' ]
l = len( let )

#Function for calculating modified z-score:
def outliers_modified_z_score( data ):

    threshold = 3.5

    median_y = np.median( data )
    median_absolute_deviation_y = np.median( [ np.abs( y - median_y ) for y in data ] )
    modified_z_scores = [ 0.6745 * ( y - median_y ) / median_absolute_deviation_y for y in data ]

    return np.where( np.abs( modified_z_scores ) > threshold )

for z in range( 5 ):

	print opts[z]

	#Load indicator files
	f1 = netcdf.netcdf_file( files[z], 'r' )

	if z < 2: #CMIP
		m = 26
		inds = np.zeros( ( l, m ) )
	elif z > 1 and z < 4: #AMIP
		m = 24
		inds = np.zeros( ( l, m ) )
	else: #AQUA
		m = 8
		inds = np.zeros( ( l, m ) )

	#Load indicators
	for i in range( l ):
		inds[i] = f1.variables[ let[i] ][:]

	#Search for outlers
	for i in range( l ):
		if len(outliers_modified_z_score( inds[i] )) > 0:
			print let[i], outliers_modified_z_score( inds[i] ) 
			





				
