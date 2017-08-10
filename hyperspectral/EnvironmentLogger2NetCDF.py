#!/usr/bin/env python

import os
import sys
import bisect
import os.path
import json
import numpy as np
from datetime import date, datetime, timedelta
from netCDF4 import Dataset



def main(elFile, ncFile):

    # read all of json in one go
    fp = open(elFile, "r")
    lines = fp.readlines()
    slines = "".join(lines)
    js = json.loads(slines)
    fp.close()

    #time stamp stuff - convert time stamp into CF Convention one
    # e.g  "timestamp": "2017.07.31-12:10:56" to "seconds since 2017-07-31 12:10:56"
    firstTime=str()
    firstTime=js["environment_sensor_readings"][0]["timestamp"]
    firstTime=firstTime.replace("."," ").replace("-"," ")
    ArrayTime=firstTime.split(" ")

    if len(ArrayTime) < 4:
        print("Cannot Parse first timestamp {0}".format(js["environment_sensor_readings"][0]["timestamp"]))
        exit(1)

    sinceTime="seconds since {0}-{1}-{2} {3}".format( ArrayTime[0],ArrayTime[1],ArrayTime[2],ArrayTime[3] )



    # create dims and coordinate vars
    output_handle=Dataset(ncFile, "w", format="NETCDF4")

    # nb time is unlimited
    output_handle.createDimension("time", 0)
    nTime=output_handle.createVariable("time", "f8", ("time",))
    setattr(nTime, "units", sinceTime)

    nSunDirection=output_handle.createVariable("sunDirection", "f8", ("time",))
    setattr(nSunDirection, "units", "degrees")

    # get wavelength  from first record only as it is repeated
    wavelength=js["environment_sensor_readings"][0]["spectrometer"]["wavelength"]
    output_handle.createDimension("wvl_lgr", len(wavelength))
    output_handle.createVariable("wvl_lgr", "f8", ("wvl_lgr",))
    # convert units to meteers
    output_handle.variables["wvl_lgr"][:] = [fd * 1e-9 for fd in wavelength]
    setattr(output_handle.variables["wvl_lgr"], "units", "meter")


    nSpectrum=output_handle.createVariable("weighted_average_downwelling_irradiance", "f8", ("time","wvl_lgr",))
    setattr(nSpectrum, "units", "watt meter-2 meter-1")

    sz = len(js["environment_sensor_readings"])

    itime=0.0 # assume that time increaes in 5 second increments across records

    for idx in range(sz):

        # output_handle.variables["time"][idx]=itime
        nTime[idx] = itime
        # sun Direction
        nSunDirection[idx]=js["environment_sensor_readings"][idx]["weather_station"]["sunDirection"]["value"]

        #
        spectrum=js["environment_sensor_readings"][idx]["spectrometer"]["spectrum"]
        nSpectrum[idx,:] = spectrum

        # assume that time increaes in 5 second increments across records
        itime+=5.0

    output_handle.close()

if __name__ == "__main__":

    argc=len(sys.argv)
    if argc == 1  :
        print("Please specify a json input file")
        exit(1)
    elif argc == 2 :
        ncFile=sys.argv[1].replace(".json",".nc")
        main(sys.argv[1],ncFile)
    elif argc >2:
        main(sys.argv[1], sys.argv[2]);




