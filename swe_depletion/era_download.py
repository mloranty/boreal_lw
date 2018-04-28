#!/usr/bin/env python

from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
   
# script to download era-interim climare reanlaysis data
# you must haver an account with ECMWFDataServer
# for info on era-interim see https://www.ecmwf.int/sites/default/files/elibrary/2011/8174-era-interim-archive-version-20.pdf
# for example python scripts see https://software.ecmwf.int/wiki/display/WEBAPI/Python+ERA-interim+examples
   
server.retrieve({
    # Specify the ERA-Interim data archive. Don't change.
	'class'     : "ei",
	'stream'    : "oper",
	# analysis (an) vars, rather than forecast
	'type'      : "an",
	# surface variables (as opposed to atm pressure levels)
    'levtype'   : "sfc",
	# select parameters, for codes see http://apps.ecmwf.int/codes/grib/param-db
	# this is 2m air temp
    'param'     : "167.128",
    'dataset'   : "interim",
	# timestep info - see also https://www.ecmwf.int/en/faq/what-are-steps-surface-daily-fields-era-interim
	# and see https://software.ecmwf.int/wiki/pages/viewpage.action?pageId=56658233
    'step'      : "0",
    'time'      : "00/06/12/18",
    'date'      : "2015-02-01/to/2015-06-30",
	#specify resolution in degrees lat/lon
	'grid'      : "0.5/0.5",
	# optionally restrict area (in N/W/S/E).
    'area'      : "90/-180/50/180",
	# file format and output
    'format'    : "netcdf",
    'target'    : "interim_2015-02-01to2015-06-30_00061218.nc"
 })
 
 
 

 

