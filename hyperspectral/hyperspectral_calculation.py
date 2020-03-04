# -*- coding: utf-8 -*-
import numpy as np
import os
import json
from math import *
from datetime import date, datetime, timedelta
from decimal import *

from terrautils.spatial import scanalyzer_to_latlon


LATITUDE_TO_METER = 1 / (30.87 * 3600)
LONGITUDE_TO_METER  = 1 / (25.906 * 3600) # varies, but has been corrected based on the actual location of the field


# --- TESTING FUNCTIONS - DON'T HARDCODE FOR PRODUCTION ---
def get_VNIR_fixed_metadata(date):
    """Normally, this is included in the cleaned metadata that terrautils will generate, using the Github fixed
        metadata repo as a source.

        I'm copying this code to this test file for convenience so we can test using the raw metadata file without
        cleaning it first. Only fields that this code might use is included.

        Numbers from: https://github.com/terraref/sensor-metadata/blob/master/sensors/VNIR/sensor_fixed_metadata.json
        Cleaning code: https://github.com/terraref/terrautils/blob/master/terrautils/lemnatec.py#L169
    """

    # New camera
    if date >= "2018-08-18":
        md = {
            "location_in_camera_box_m": {
                "x": "0.867",
                "y": "1.663",
                "z": "0.635"
            },
            "field_of_view_degrees": {
                "y": "21.0"
            }
        }
    # Old camera
    else:
        md = {
            "location_in_camera_box_m": {
                "x": "0.877",
                "y": "2.325",
                "z": "0.635"
            },
            "field_of_view_degrees": {
                "y": "21.0"
            }
        }

    return md

def get_SWIR_fixed_metadata(date):
    """Normally, this is included in the cleaned metadata that terrautils will generate, using the Github fixed
        metadata repo as a source.

        I'm copying this code to this test file for convenience so we can test using the raw metadata file without
        cleaning it first. Only fields that this code might use is included.

        Numbers from: https://github.com/terraref/sensor-metadata/blob/master/sensors/SWIR/sensor_fixed_metadata.json
        Cleaning code: https://github.com/terraref/terrautils/blob/master/terrautils/lemnatec.py#L169
    """

    # New camera
    if date >= "2018-08-18":
        md = {
            "location_in_camera_box_m": {
                "x": "0.495",
                "y": "1.875",
                "z": "0.635"
            },
            "field_of_view_degrees": {
                "y": "44.5"
            }
        }
    # Old camera
    else:
        md = {
            "location_in_camera_box_m": {
                "x": "0.877",
                "y": "2.325",
                "z": "0.635"
            },
            "field_of_view_degrees": {
                "y": "44.5"
            }
        }

    return md

def generate_geojson(SE, SW, NE, NW):
    # Create a geojson file for viewing e.g. in QGIS
    return {
        "type": "Polygon",
        "coordinates": [[
            [NE[1], NE[0]],
            [NW[1], NW[0]],
            [SW[1], SW[0]],
            [SE[1], SE[0]],
            [NE[1], NE[0]] # repeat final coordinate
        ]]
    }

def test_pixel2Geographic(json_path, hdr_path, camera_type):
    # operates off raw json md file instead of cleaned metadata in pipeline
    synthetic_metadata = {}

    # ----- GET FIXED METADATA -----
    date = hdr_path.split("/")[-3]
    timestamp = hdr_path.split("/")[-2]
    if camera_type == "SWIR":
        synthetic_metadata['sensor_fixed_metadata'] = get_SWIR_fixed_metadata(date)
    else:
        synthetic_metadata['sensor_fixed_metadata'] = get_VNIR_fixed_metadata(date)

    # ----- GET VARIABLE METADATA -----
    with open(json_path) as json_meta:
        meta = json.loads(json_meta.read())["lemnatec_measurement_metadata"]
    x_gantry_pos, y_gantry_pos = None, None
    # For scans east-to-west, the captured position is in the southeast.
    # For scans west-to-east, the captured position is the southwest.
    for pos_field in ["position x [m]", "Position x [m]"]:
        if pos_field in meta["gantry_system_variable_metadata"]:
            x_gantry_pos = float(meta["gantry_system_variable_metadata"][pos_field])
            y_gantry_pos = float(meta["gantry_system_variable_metadata"][pos_field.replace("x", "y")])
            break

    synthetic_metadata['gantry_variable_metadata'] = {
        'position_m': {
            'x': x_gantry_pos,
            'y': y_gantry_pos
        },
        'scan_direction_is_positive': str(meta["gantry_system_variable_metadata"]["scanDirectionIsPositive"]),
        'scan_speed_m/s': float(meta["gantry_system_variable_metadata"]["scanSpeedInMPerS [m/s]"])
    }
    synthetic_metadata['sensor_variable_metadata'] = {
        'frame_period': int(meta["sensor_variable_metadata"]["current setting frameperiod"])
    }

    out = pixel2Geographic(synthetic_metadata, hdr_path, camera_type)

    out_path = "/home/extractor/hs_calib/VNIR/2019-07-25/geojson/%s.geojson" % timestamp
    with open(out_path, 'w') as shpfile:
        json.dump(out['bbox_geojson'], shpfile)

    return out

# ----------------- END TESTING FUNCTIONS -----------------

def solar_zenith_angle(time_date):

    latitude = 33 + 4.47 / 60

    days_offset       = time_date - datetime(year=time_date.year,month=1,day=1) + \
                        timedelta(days=1)
    numerical_cal_day = Decimal(days_offset.days) + Decimal(days_offset.seconds) / \
                        Decimal(86340)

    theta = Decimal(2)*Decimal(pi)*numerical_cal_day/Decimal(365)

    eccentricity_fac = Decimal(1.000110)+ \
                       Decimal(cos(theta))*Decimal(0.034221)+ \
                       Decimal(sin(theta))*Decimal(0.001280)+ \
                       Decimal(cos(Decimal(2)*theta))*Decimal(0.000719)+ \
                       Decimal(sin(Decimal(2)*theta))*Decimal(0.000077)

    solar_decline    = Decimal(0.006918) - Decimal(0.399912)*Decimal(cos(theta))+ \
                       Decimal(0.070257)*Decimal(sin(theta))- \
                       Decimal(0.006758)*Decimal(cos(Decimal(2)*theta))+ \
                       Decimal(0.000907)*Decimal(sin(Decimal(2)*theta))- \
                       Decimal(0.002697)*Decimal(cos(Decimal(3)*theta))+ \
                       Decimal(0.001480)*Decimal(sin(Decimal(3)*theta))

    sin_latitude     = Decimal(sin(radians(latitude)))
    sin_delta        = Decimal(sin(solar_decline))
    cos_latitude     = Decimal(cos(radians(latitude)))
    cos_delta        = Decimal(cos(solar_decline))

    cphase = Decimal(cos(Decimal(2*pi)*numerical_cal_day))
    cos_solar_zen_ang = sin_latitude*sin_delta-cos_latitude*cos_delta*cphase

    return float(Decimal(acos(cos_solar_zen_ang)/pi)*Decimal(180))

def pixel2Geographic(metadata, hdr_path, camera_type):

    # ----- PREPARE METADATA -----
    # TODO: store this in fixed metadata
    if camera_type == "SWIR":
        x_pixel_size = 1.930615052e-3
    else:
        x_pixel_size = 1.025e-3
    y_pixel_size = 0.98526434004512529576754637665e-3

    if 'scan_speed_m/s' in metadata['gantry_variable_metadata']:
        scan_speed = metadata['gantry_variable_metadata']['scan_speed_m/s']
        framepd = metadata['sensor_variable_metadata']['frame_period']
        downsampled = (scan_speed == 0.04 and framepd > 25)
    else:
        downsampled = False

    x_camera_pos = float(metadata['sensor_fixed_metadata']["location_in_camera_box_m"]["x"])
    y_camera_pos = float(metadata['sensor_fixed_metadata']["location_in_camera_box_m"]["y"])
    x_gantry_pos = float(metadata['gantry_variable_metadata']["position_m"]["x"])
    y_gantry_pos = float(metadata['gantry_variable_metadata']["position_m"]["y"])
    scan_dir = str(metadata['gantry_variable_metadata']["scan_direction_is_positive"])

    # ----- GET HDR METADATA -----
    with open(hdr_path) as hdr_data:
        overall = hdr_data.readlines()
    for members in overall:
        if "samples" in members:
            x_pixel_num = int(members.split("=")[-1].strip("\n"))
        if "lines" in members:
            y_pixel_num = int(members.split("=")[-1].strip("\n"))

    # ----- GET ORIGIN POINT -----
    REFERENCE_POINT = scanalyzer_to_latlon(0,0) # SE corner of the field latlon coordinates
    origin_x = REFERENCE_POINT[0]
    origin_y = REFERENCE_POINT[1]

    # ----- CALCULATE POSITIONS -----
    x_absolute_pos = x_gantry_pos + x_camera_pos
    y_absolute_pos = y_gantry_pos + y_camera_pos

    x_final_result = []
    x_half_px = x_pixel_num/2
    for x in range(x_pixel_num):
        if x <= x_half_px:
            x_final_result_value = x_absolute_pos + ((-x+x_half_px) * -x_pixel_size)
        else:
            x_final_result_value = x_absolute_pos + ((x-x_half_px) * x_pixel_size)
        x_final_result.append(x_final_result_value)
    x_final_result = np.array(x_final_result)

    if scan_dir == "True":
        if not downsampled:
            y_final_result = np.array([y * y_pixel_size for y in range(y_pixel_num)]) + y_absolute_pos
        else:
            y_final_result = np.array([y * 2 * y_pixel_size for y in range(y_pixel_num)]) + y_absolute_pos
    else:
        if not downsampled:
            y_final_result = np.array([-y * y_pixel_size for y in range(y_pixel_num)]) + y_absolute_pos
        else:
            y_final_result = np.array([-y * 2 * y_pixel_size for y in range(y_pixel_num)]) + y_absolute_pos

    # Determine 4 corners of bounding box
    SE = (x_final_result[0]  * LATITUDE_TO_METER + origin_x, -(y_final_result[0])  * LONGITUDE_TO_METER + origin_y)
    SW = (x_final_result[0]  * LATITUDE_TO_METER + origin_x, -(y_final_result[-1]) * LONGITUDE_TO_METER + origin_y)
    NE = (x_final_result[-1] * LATITUDE_TO_METER + origin_x, -(y_final_result[0])  * LONGITUDE_TO_METER + origin_y)
    NW = (x_final_result[-1] * LATITUDE_TO_METER + origin_x, -(y_final_result[-1]) * LONGITUDE_TO_METER + origin_y)

    j_bbox = generate_geojson(SE, SW, NE, NW)

    lat_final_result = np.array([x * LATITUDE_TO_METER  for x in x_final_result]) + origin_x
    lon_final_result = np.array([-y * LONGITUDE_TO_METER for y in y_final_result]) + origin_y
    bounding_box = (str(SE).strip("()"), str(SW).strip("()"), str(NE).strip("()"), str(NW).strip("()"))

    return {"x_coordinates": x_final_result,
            "y_coordinates": y_final_result,
            "x_pixel_size" : x_pixel_size,
            "y_pixel_size" : y_pixel_size,
            "x_reference"  : origin_x,
            "y_reference"  : origin_y,
            "latitudes"    : lat_final_result,
            "longitudes"   : lon_final_result,
            "bounding_box" : bounding_box,
            "bbox_geojson" : j_bbox}
