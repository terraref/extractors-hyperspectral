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
    '''
    Calculate the solar zenith angle based on the datetime object given
    '''
    
    latitude  = Decimal(33)  + Decimal(4.47)   / Decimal(60)
    longitude = Decimal(-111)- Decimal(58.485) / Decimal(60)
    
    getcontext().prec = 8
    julian_date = _julian_date(time_date)
    actural_day = Decimal(julian_date)
    
    ### Sun & Earth Data ###
    long_perihelion = Decimal(282.9404) + Decimal(4.70935e-5) * actural_day
    eccentricity    = Decimal(0.016709) - Decimal(1.151e-9)   * actural_day
    anomaly_degree  = (Decimal(356.047) + Decimal(0.9856002585)*actural_day) % Decimal(360)
    
    sun_mean_long   = long_perihelion + anomaly_degree
    sun_obliquity   = Decimal(23.4393)- Decimal(3.563e-7) * actural_day
    auxiliary_angle = anomaly_degree + (Decimal(180)/Decimal(pi)) * eccentricity *\
                      Decimal(sin(anomaly_degree*Decimal(pi)/Decimal(180))) *\
                      (Decimal(1)+eccentricity*Decimal(cos(anomaly_degree*Decimal(pi)/Decimal(180))))
            
    x_rectan_coord  = Decimal(cos(auxiliary_angle*Decimal(pi/180))) - eccentricity
    y_rectan_coord  = Decimal(sin(auxiliary_angle*Decimal(pi/180))) * Decimal(sqrt(float(Decimal(1)-\
                      eccentricity**Decimal(2))))    

    distance        = Decimal(sqrt(float(x_rectan_coord**2+y_rectan_coord**2)))
    true_anomaly    = Decimal(atan2(y_rectan_coord, x_rectan_coord))*Decimal(180/pi)
    sun_longitude   = true_anomaly + long_perihelion
    
    x_ecliptic      = distance*Decimal(cos(sun_longitude*Decimal(pi/180)))
    y_ecliptic      = distance*Decimal(sin(sun_longitude*Decimal(pi/180)))
    z_ecliptic      = Decimal(0.0)
    
    x_equitorial    = x_ecliptic
    y_equitorial    = y_ecliptic*Decimal(cos(sun_obliquity*Decimal(pi/180)))+\
                      z_ecliptic*Decimal(sin(sun_obliquity*Decimal(pi/180)))
    z_equitorial    = y_ecliptic*Decimal(sin(Decimal(23.4406)*Decimal(pi/180)))+\
                      z_ecliptic*Decimal(cos(sun_obliquity*Decimal(pi/180)))
    
    distance        = (x_equitorial**2+y_equitorial**2+z_equitorial**2)**(1/2)
    delta           = Decimal(asin(z_equitorial/distance))*Decimal(180/pi)
    right_ascension = Decimal(atan2(y_equitorial, x_equitorial))*Decimal(180/pi)
    
    uth             = Decimal(time_date.hour)+Decimal(time_date.minute)/Decimal(60)+\
                      Decimal(time_date.second)/Decimal(3600)
    gmst            = ((sun_longitude+Decimal(180))%Decimal(360))/Decimal(15)
    loc_sideral_t   = gmst + uth + longitude / Decimal(15)

    ### Right Ascension to Hour Angle ###
    hour_angle      = loc_sideral_t*Decimal(15) - right_ascension
    
    actural_x       = Decimal(cos(hour_angle*Decimal(pi/180)))*Decimal(cos(delta*Decimal(pi/180)))
    actural_z       = Decimal(sin(delta*Decimal(pi/180)))
    
    x_rotation      = actural_x*Decimal(cos((Decimal(90)-latitude)*Decimal(pi)/Decimal(180)))-\
                      actural_z*Decimal(sin((Decimal(90)-latitude)*Decimal(pi)/Decimal(180)))
    z_rotation      = actural_x*Decimal(sin((Decimal(90)-latitude)*Decimal(pi/180)))+\
                      actural_z*Decimal(cos((Decimal(90)-latitude)*Decimal(pi/180)))

    ### Zenith Angle = 90 - Elevation Angle ###  
    return 90 - abs(Decimal(asin(z_rotation))*Decimal(180/pi))


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