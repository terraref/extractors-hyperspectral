#!/usr/bin/env python

import os
import sys
import os.path
from datetime import datetime, timedelta

'''
The date should be in format
month/date/year hour:minute:second (to be compatible with the metadata in the hyperspectral workflow)
The path to the nearest two environmental logger records (can be JSON or netCDF file) will be printed out to the STDOUT

Sample Command-line (on ROGER):
python hyperspectral_file_finder.py "09/22/2016 13:53:44" /projects/arpae/terraref/sites/ua-mac/raw_data/EnvironmentLogger

Generally the commands look like:
python hyperspectral_file_finder.py <"month/date/year hour:minute:second"> <root directory of the environmental logger>

# Make sure that you quote the time parameter.
'''

def date_parser(date_object, time_shift_unit, time_shift):
    current_time = date_object
    if time_shift_unit == "hours":
        next_period = date_object + timedelta(hours=time_shift)
    elif time_shift_unit == "minutes":
        next_period = date_object + timedelta(minutes=time_shift)

    return {"current_time_year" : str(current_time.year),
            "current_time_month": format(current_time.month, "02"),
            "current_time_day"  : format(current_time.day, "02"),
            "current_time_hour" : format(current_time.hour,"02"),
            "next_period_year"  : str(next_period.year),
            "next_period_month" : format(next_period.month, "02"),
            "next_period_day"   : format(next_period.day, "02"),
            "next_period_hour"  : format(next_period.hour,"02")}


def analyze_time_interval(target, date_a, date_b):

def file_finder(date, file_path):
    '''
    Will help find out the most close environmental logger files
    '''
    if not os.path.exists(file_path):
        raise OSError("The path is invalid")

    date_object = datetime.strptime(date, "%m/%d/%Y %H:%M:%S")
    current_date_directory = "-".join((str(date_object.year), format(date_object.month, "02"), format(current_time.day, "02")))

    if len(os.listdir(os.path.join(file_path, current_date_directory))) > 24:
        date_list = date_parser(date_object,"minute",2)
    current_date_directory = "-".join((date_list["current_time_year"], date_list["current_time_month"], date_list["current_time_day"]))
    for files in os.listdir(os.path.join(file_path, current_date_directory)):
        if files.startswith(current_date_directory+"_"+date_list["current_time_hour"]): #[0][-29:]
            print os.path.join(file_path, current_date_directory, files)
            break

    next_date_directory = "-".join((date_list["next_period_year"], date_list["next_period_month"], date_list["next_period_day"]))
    for files in os.listdir(os.path.join(file_path, next_date_directory)):
        if files.startswith(next_date_directory+"_"+date_list["next_period_hour"]): #[0][-29:]
            print os.path.join(file_path, next_date_directory, files)
            break


if __name__ == "__main__":
    file_finder(sys.argv[1], sys.argv[2])