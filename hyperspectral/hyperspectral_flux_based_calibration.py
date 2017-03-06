#!/usr/bin/env python

import os
import sys
import bisect
import os.path
import numpy as np
from datetime import date, datetime, timedelta
from netCDF4 import Dataset

'''
The date should be in format
month/date/year hour:minute:second (to be compatible with the metadata in the hyperspectral workflow)
The path to the nearest two environmental logger records (can be JSON or netCDF file) will be printed out to the STDOUT

Sample Command-line (on ROGER):
python hyperspectral_flux_based_calibration.py "09/22/2016 13:53:44" /projects/arpae/terraref/sites/ua-mac/raw_data/EnvironmentLogger

Generally the commands look like:
python hyperspectral_flux_based_calibration.py <"month/date/year hour:minute:second"> <root directory of the environmental logger>

# Make sure that you quote the time parameter.

For the environmental logger records after 04/26/2016, there would be 24 files per day (1 file per hour, 5 seconds per record), the loggers before this
time point would be 720 files per day (1 file per 2 minutes, 7 seconds per record)
'''
class CalibrationError(Exception):
    pass


def solar_zenith_angle(datetime_object):
    raise NotImplementedError("Haven't finished yet'")


def weighted_avg_calculator(target_time, time_a, time_b, flx_a, flx_b):
    '''
    Calculate the weighted average for the two environmental loggers around the hyperspectral time point
    '''
    target_time_absolute = translate_time(target_time)

    weight_a = abs(time_a-target_time_absolute)
    weight_b = abs(time_b-target_time_absolute)

    #### Weighted average = sum(weight[i]*sample[i]) / sum(weights)

    return (np.array(flx_a)*weight_a + np.array(flx_b)*weight_b)/(weight_a+weight_b)


def translate_time(datetime_object):
    '''
    translate the datetime object in terms of days offset to the unix base time (01/01/1970)
    '''
    datetime_object_unpack = datetime_object.timetuple()
    unix_basetime          = date(year=1970, month=1, day=1)

    time_split  = date(year=datetime_object_unpack.tm_year, 
                       month=datetime_object_unpack.tm_mon, 
                       day=datetime_object_unpack.tm_mday) - unix_basetime

    return (time_split.total_seconds() + datetime_object_unpack.tm_hour * 3600.0 + datetime_object_unpack.tm_min * 60.0 + datetime_object_unpack.tm_sec) / (3600.0 * 24.0)


def downwelling_irradiance_extractor(netCDF_handles, target_time):
    '''
    Extract the downwelling_irradiance_spectrum for the certain time points in a netCDF file.
    '''
    if len(netCDF_handles) == 1:
        with Dataset(netCDF_handles[0], "r", format='NETCDF4') as netCDF_handler:
            numerical_time = translate_time(target_time)
            nearest_indices = bisect.bisect_left(netCDF_handler.variables["time"], numerical_time) # right time point

            return netCDF_handler.variables["time"][nearest_indices],\
                   netCDF_handler.variables["time"][nearest_indices-1],\
                   netCDF_handler.variables["flx_spc_dwn"][nearest_indices],\
                   netCDF_handler.variables["flx_spc_dwn"][nearest_indices-1]
    else:
        with Dataset(netCDF_handles[0], "r", format='NETCDF4') as netCDF_handler_a,\
             Dataset(netCDF_handles[1], "r", format='NETCDF4') as netCDF_handler_b:

            return netCDF_handler_a.variables["flx_spc_dwn"][-1], netCDF_handler_b.variables["flx_spc_dwn"][0]

def wavelength_extractor(netCDF_handles):
    '''
    Extract the wavelengths from the environmental loggers
    '''
    if len(netCDF_handles) == 1:
        with Dataset(netCDF_handles[0], "r", format='NETCDF4') as netCDF_handler:
            wavelength_list = list()
            wavelength_list[:] = netCDF_handler.variables["wvl_lgr"] # Strange, but we cannot return a variable after the dataset was closed

            return wavelength_list


def file_locator(date, file_path):
    '''
    Will help find out the most close environmental logger file
    '''
    if not os.path.exists(file_path):
        raise CalibrationError("The environmental logger root path is invalid, HINT: check your directory input")

    current_date_directory = os.path.join(file_path, "-".join((str(date.year), format(date.month, "02"), format(date.day, "02"))))

    if not os.path.exists(current_date_directory):
        raise CalibrationError("Environmental Logger does not have the record in this certain day")

    date_cutpoint = datetime(year=2016,month=4,day=26) # Since then the environmental logger generate 24 files per day

    if date > date_cutpoint:

        temp_file_header = "-".join((str(date.year), format(date.month, "02"), format(date.day, "02"))) + str(date.hour)
        directory = os.listdir(current_date_directory)
        file_time_list = [datetime.strptime(files[:19], "%Y-%m-%d_%H-%M-%S") for files in directory]

        intended_file_location = bisect.bisect_left(file_time_list, date)

        if intended_file_location > 0:
            return [os.path.join(current_date_directory, directory[intended_file_location-1])]


def main(date_string, EL_root_directory):
    date_object = datetime.strptime(date_string, "%m/%d/%Y %H:%M:%S") # parse the time
    file_list   = file_locator(date_object, EL_root_directory) # locate the file

    flx_spc_dwn = downwelling_irradiance_extractor(file_list, date_object) # extarct the downwelling irradiance data
    wavelength  = wavelength_extractor(file_list)

    ########## Write into the output file ##########
    with Dataset("calibration_el_realtime.nc","w",format="NETCDF4") as output_handle:
        output_handle.createDimension("wvl", len(wavelength))
        output_handle.createVariable("weighted_average_downwelling_irradiance", "f8", ("wvl",))
        output_handle.createVariable("wvl_lgr", "f8", ("wvl",))
        result = weighted_avg_calculator(date_object, *flx_spc_dwn)

        output_handle.variables["wvl_lgr"][:] = wavelength
        output_handle.variables["weighted_average_downwelling_irradiance"][:] = result

        setattr(output_handle.variables["weighted_average_downwelling_irradiance"], "units", "watt meter-2 meter-1")
        setattr(output_handle.variables["wvl_lgr"], "units", "meter")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])