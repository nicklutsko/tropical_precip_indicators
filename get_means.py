"""

Script to get means of indicators for CMIP, AMIP and AQUA simulations

Copyright Nick Lutsko and Max Popp 2017

"""
from scipy.io import netcdf
import numpy as np

print "Calculate means"

#Different sets of simulations
opts = [ "CMIP Historical 1979 - 2005", "CMIP Historical 1861 - 2005", "AMIP 1979 - 2005", "AMIP 1979 - 2008", "aqua 1979 - 1983 " ]

#File names
files = [ 'Precipitation_indicators_SAH_historical_r1i1p1_197901-200512_ZM.nc', 'Precipitation_indicators_SAH_historical_r1i1p1_186101-200512_ZM.nc', 'Precipitation_indicators_SAH_amip_r1i1p1_197901-200512_ZM.nc', 
'Precipitation_indicators_SAH_amip_r1i1p1_197901-200812_ZM.nc',
'Precipitation_indicators_SAH_aqua_r1i1p1_197901-198312_ZM.nc' ]

#Indicators
let = [ 'E', 'D', 'A', 'H', 'R', 'I', 'S', 'L', 'W' ]

obs = [ 4.02, 1.88, 1.81, 4.20, -0.15, 6.25, 10., 17.5, 45. ] #GPCP values
obs2 = [ 4.95, 2.62, 2.26, 5.51, 0.17, 6.75, 11.25, 18.00, 48.75 ] #ERA-Interim values

for z in range( 5 ):
	print opts[z]

	#Load file
	f1 = netcdf.netcdf_file( files[z], 'r' )

	if z < 2:
		#CMIP simulations
		inds = np.zeros( ( 9, 26 ) )
	elif z > 1 and z < 4:
		#AMIP simulations
		inds = np.zeros( ( 9, 24 ) )
	else:
		#AQUA simulations
		inds = np.zeros( ( 9, 8 ) )

	#Load indicators
	for i in range( 9 ):
		inds[i] = f1.variables[ let[i] ][:]

	
	for i in range( 9 ):
		a = 0
		#Print means, alerting if means are 2 std from obs
		if np.mean( inds[i] ) + 2. * np.std( inds[i] ) <  obs[i]:
			print let[i], "GPCP", np.mean( inds[i] ), obs[i] 	
			a += 1
		if np.mean( inds[i] ) - 2. * np.std( inds[i] ) >  obs[i] :
			print let[i], "GPCP", np.mean( inds[i] ), obs[i]
			a += 1
		if  np.mean( inds[i] )  + 2. * np.std( inds[i] ) < obs2[i] :
			print let[i], "ERA", np.mean( inds[i] ),  obs2[i]
			a += 1
		if np.mean( inds[i] )  - 2. * np.std( inds[i] ) > obs2[i] :
			print let[i], "ERA", np.mean( inds[i] ),  obs2[i]
			a += 1
		if a == 0:
			print let[i], np.mean( inds[i] )

