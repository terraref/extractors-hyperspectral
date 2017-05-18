# -*- coding: utf-8 -*-

import numpy as np
import sys
import json
from math import *
from datetime import date, datetime, timedelta
from decimal import *

# from Dr. LeBauer, Github thread: terraref/referece-data #32
CAMERA_POSITION = np.array([1.9, 0.855, 0.635])

# from Dr. LeBauer, Github thread: terraref/referece-data #32
CAMERA_FOCAL_LENGTH = 24e-3 # the focal length for SWIR camera. unit:[m]

# from Dr. LeBauer, Github thread: terraref/referece-data #32
PIXEL_PITCH = 25e-6 #[m]

REFERENCE_POINT = 33 + 4.47 / 60, -111 - 58.485 / 60 # from https://github.com/terraref/reference-data/issues/32

LONGITUDE_TO_METER = 1 / (30.87 * 3600)
LATITUDE_TO_METER  = 1/ (25.906 * 3600) #varies, but has been corrected based on the actural location of the field

GOOGLE_MAP_TEMPLATE = "https://maps.googleapis.com/maps/api/staticmap?size=1280x720&zoom=17&path=color:0x0000005|weight:5|fillcolor:0xFFFF0033|{pointA}|{pointB}|{pointC}|{pointD}"


def _julian_date(time_date):
    '''
    Calculate the Julian Date based on the input datetime object (should
    be in Gregorian).
    
    Private in this module

    Had already checked the output against the result from
    United States Naval Observatory, Astronomical App Dept.
    
    Detail: http://aa.usno.navy.mil/data/docs/JulianDate.php
    '''
    a = floor((14-time_date.month)/12)
    
    if time_date.month in (1, 2):
        assert a == 1
    else:
        assert a == 0
        
    years  = time_date.year + 4800 - a
    months = time_date.month + 12 * a - 3
    
    if time_date.month == 3:
        assert months == 0
    elif time_date.month == 2:
        assert months == 11
    
    getcontext().prec = 8
    return time_date.day  + floor((153*months+2)/5) + 365*years +\
           floor(years/4) - floor(years/100) + floor(years/400) - 32045 +\
           float(Decimal(time_date.hour-12)/Decimal(24) + Decimal(time_date.minute)/Decimal(1440)+\
           Decimal(time_date.second)/Decimal(86400))
        
        
def solar_zenith_angle(time_date):

    latitude = 33 + 4.47 / 60
    
    days_offset       = time_date - datetime(year=time_date.year,month=1,day=1) +\
                        timedelta(days=1)
    numerical_cal_day = Decimal(days_offset.days) + Decimal(days_offset.seconds) /\
                        Decimal(86340)
    
    theta = Decimal(2)*Decimal(pi)*numerical_cal_day/Decimal(365)
    
    eccentricity_fac = Decimal(1.000110)+\
                       Decimal(cos(theta))*Decimal(0.034221)+\
                       Decimal(sin(theta))*Decimal(0.001280)+\
                       Decimal(cos(Decimal(2)*theta))*Decimal(0.000719)+\
                       Decimal(sin(Decimal(2)*theta))*Decimal(0.000077)
    
    solar_decline    = Decimal(0.006918) - Decimal(0.399912)*Decimal(cos(theta))+\
                                           Decimal(0.070257)*Decimal(sin(theta))-\
                                           Decimal(0.006758)*Decimal(cos(Decimal(2)*theta))+\
                                           Decimal(0.000907)*Decimal(sin(Decimal(2)*theta))-\
                                           Decimal(0.002697)*Decimal(cos(Decimal(3)*theta))+\
                                           Decimal(0.001480)*Decimal(sin(Decimal(3)*theta))
                        
    sin_latitude     = Decimal(sin(radians(latitude)))
    sin_delta        = Decimal(sin(solar_decline))
    cos_latitude     = Decimal(cos(radians(latitude)))
    cos_delta        = Decimal(cos(solar_decline))
    
    cphase = Decimal(cos(Decimal(2*pi)*numerical_cal_day))
    cos_solar_zen_ang = sin_latitude*sin_delta-cos_latitude*cos_delta*cphase
    
    return float(Decimal(acos(cos_solar_zen_ang)/pi)*Decimal(180))


def pixel2Geographic(jsonFileLocation, headerFileLocation, cameraOption, downsampled=False):

    ######################### Load necessary data #########################
    with open(jsonFileLocation) as fileHandler:
        master = json.loads(fileHandler.read())["lemnatec_measurement_metadata"]

        if "position x [m]" in master["gantry_system_variable_metadata"]:
            x_gantry_pos = float(master["gantry_system_variable_metadata"]["position x [m]"])
            y_gantry_pos = float(master["gantry_system_variable_metadata"]["position y [m]"])

        elif "Position x [m]" in master["gantry_system_variable_metadata"]:
            x_gantry_pos = float(master["gantry_system_variable_metadata"]["Position x [m]"])
            y_gantry_pos = float(master["gantry_system_variable_metadata"]["Position y [m]"])

        else: # We notice that there are cases that no position data available
            return {"x_coordinates": None,
                    "y_coordinates": None,
                    "latitudes"    : None,
                    "longitudes"   : None,
                    "bounding_box" : None,
                    "Google_Map"   : None}

        x_camera_pos = 1.9 # From https://github.com/terraref/reference-data/issues/32
        y_camera_pos = 0.855

        if cameraOption == "SWIR":
            x_pixel_size = 1.930615052e-3
        else:
            x_pixel_size = 1.025e-3

        y_pixel_size = 0.98526434004512529576754637665e-3

        with open(headerFileLocation) as fileHandler:
            overall = fileHandler.readlines()

            for members in overall:
                if "samples" in members:
                    x_pixel_num = int(members.split("=")[-1].strip("\n"))
                elif "lines" in members:
                    y_pixel_num = int(members.split("=")[-1].strip("\n"))


        ######################### Do calculation #########################

        x_absolute_pos = x_gantry_pos + x_camera_pos
        y_absolute_pos = y_gantry_pos + y_camera_pos

        x_final_result = np.array([x * x_pixel_size for x in range(x_pixel_num)]) + x_absolute_pos

        if not downsampled:
            y_final_result = np.array([y * y_pixel_size for y in range(y_pixel_num)]) + y_absolute_pos
        else:
            y_final_result = np.array([y * 2 * y_pixel_size for y in range(y_pixel_num)]) + y_absolute_pos

        ########### Sample result: x -> 0.377 [m], y -> 0.267 [m] ###########

        SE = x_final_result[-1] * LONGITUDE_TO_METER + REFERENCE_POINT[0], y_final_result[-1] * LATITUDE_TO_METER + REFERENCE_POINT[1]
        SW = x_final_result[0] * LONGITUDE_TO_METER  + REFERENCE_POINT[0], y_final_result[-1] * LATITUDE_TO_METER + REFERENCE_POINT[1]
        NE = x_final_result[-1] * LONGITUDE_TO_METER + REFERENCE_POINT[0], y_final_result[0]  * LATITUDE_TO_METER + REFERENCE_POINT[1]
        NW = x_final_result[0] * LONGITUDE_TO_METER + REFERENCE_POINT[0] , y_final_result[0]  * LATITUDE_TO_METER + REFERENCE_POINT[1]

        bounding_box = (str(SE).strip("()"), str(SW).strip("()"), str(NE).strip("()"), str(NW).strip("()"))
        bounding_box_mapview = GOOGLE_MAP_TEMPLATE.format(pointA=bounding_box[0],
                                                          pointB=bounding_box[1],
                                                          pointC=bounding_box[2],
                                                          pointD=bounding_box[3])

        lat_final_result = np.array([x * LATITUDE_TO_METER for x in x_final_result])  + REFERENCE_POINT[0]
        lon_final_result = np.array([-y * LONGITUDE_TO_METER for y in y_final_result]) + REFERENCE_POINT[1]

        return {"x_coordinates": x_final_result,
                "y_coordinates": y_final_result,
                "latitudes"    : lat_final_result,
                "longitudes"   : lon_final_result,
                "bounding_box" : bounding_box,
                "Google_Map"   : bounding_box_mapview}