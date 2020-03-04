"""
*** Hyperspectral Image Calibration ***
This module will process the raw data file and export ENVI files with variables
stored its hdr file. (This is an initial version, more improvements are coming up. )

input:   1) raw VNIR data;
         2) spectral data from downwelling irradiance sensor

output: the reflectance image

@author: SLU Remote Sensing Lab
"""

import json
import numpy as np
import os
import spectral.io.envi as envi
from netCDF4 import Dataset
from PIL import Image
from datetime import date, datetime, timedelta

from hyperspectral_calculation import pixel2Geographic, test_pixel2Geographic, solar_zenith_angle


raw_root = "/home/extractor/hs_calib"
calib_root = "/home/extractor"


def create_empty_netCDF(raw_file, out_file):
    if not os.path.isdir(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))

    empty = Dataset(out_file, "w", format="NETCDF4")

    # Get values dynamically from Header file
    hdr_file = raw_file + ".hdr"
    hdr_samples = None
    hdr_lines = None
    hdr_bands = None
    with open(hdr_file, 'r') as hdr:
        for l in hdr.readlines():
            if l.startswith('samples'):
                hdr_samples = int(l.split("=")[1].strip())
            if l.startswith('lines'):
                hdr_lines = int(l.split("=")[1].strip())
            if l.startswith('bands'):
                hdr_bands = int(l.split("=")[1].strip())
            if hdr_samples is not None and hdr_lines is not None and hdr_bands is not None:
                break

    # Get number of frames
    frame_file = raw_file.replace("_raw", "_frameIndex.txt")
    with open(frame_file, 'r') as frames:
        num_frames = len(frames.readlines()[1:])

    # add dimensions
    empty.createDimension("wavelength", hdr_bands)
    empty.createDimension("x", hdr_lines)
    empty.createDimension("y", hdr_samples)
    empty.createDimension("time", num_frames)

    # raw reflectance
    # TODO: Do we need this if it is never populated?
    v = empty.createVariable("xps_img", "u2", ("wavelength", "x", "y"))
    v.long_name = "Exposure counts"
    v.meaning = "Counts on scale from 0 to 2^16-1 = 65535"
    v.units = "1"

    # calibrated reflectance
    v = empty.createVariable("rfl_img", "f8", ("wavelength", "x", "y"))
    v.long_name = "Reflectance of image"
    v.standard_name = "surface_albedo"
    v.units = "1"

    # solar zenith angle
    v = empty.createVariable("solar_zenith_angle", "f8", ("time"))
    v.long_name = "Solar Zenith Angle"
    v.units = "degree"
    v.notes = "The angle of the sun comparing to the vertical axis of the Cartesian Coordinate"
    v.acknowledgements = "Algorithm provided by Charles S. Zender, this Python implementation was translated from his original C program"

    # frametime
    v = empty.createVariable("frametime", "f8", ("time"))
    v.units = "days since 1970-01-01 00:00:00"
    v.calendar = "gregorian"
    v.notes = "date stamp per each scanline"

    # raw location
    v = empty.createVariable("x", "f8", ("x"))
    v.units = "meter"
    v.long_name = "North distance from southeast corner of field"
    v = empty.createVariable("y", "f8", ("y"))
    v.units = "meter"
    v.long_name = "West distance from southeast corner of field"
    v = empty.createVariable("latitude", "f8", ("x"))
    v.units = "degree_north"
    v.long_name = "The precise latitude value for each pixel in the picture"
    v = empty.createVariable("longitude", "f8", ("y"))
    v.units = "degree_east"
    v.long_name = "The precise longitude value for each pixel in the picture"

    # pixel sizes
    v = empty.createVariable("x_pxl_sz", "f8")
    v.long_name = "x coordinate length of a single pixel in VNIR images"
    v.units = "meters"
    v = empty.createVariable("y_pxl_sz", "f8")
    v.long_name = "y coordinate length of a single pixel in pictures captured by SWIR and VNIR camera"
    v.units = "meters"

    # gantry positional variables
    v = empty.createVariable("x_img_ne", "f8")
    v.long_name = "Northeast corner of image, north distance to reference point"
    v.units = "meters"
    v = empty.createVariable("x_img_nw", "f8")
    v.long_name = "Northwest corner of image, north distance to reference point"
    v.units = "meters"
    v = empty.createVariable("x_img_se", "f8")
    v.long_name = "Southeast corner of image, north distance to reference point"
    v.units = "meters"
    v = empty.createVariable("x_img_sw", "f8")
    v.long_name = "Southwest corner of image, north distance to reference point"
    v.units = "meters"
    v = empty.createVariable("x_reference_point", "f8")
    v.long_name = "x of the master reference point at southeast corner of field"
    v.units = "meters"
    v = empty.createVariable("y_img_ne", "f8")
    v.long_name = "Northeast corner of image, west distance to reference point"
    v.units = "meters"
    v = empty.createVariable("y_img_nw", "f8")
    v.long_name = "Northwest corner of image, west distance to reference point"
    v.units = "meters"
    v = empty.createVariable("y_img_se", "f8")
    v.long_name = "Southeast corner of image, west distance to reference point"
    v.units = "meters"
    v = empty.createVariable("y_img_sw", "f8")
    v.long_name = "Southwest corner of image, west distance to reference point"
    v.units = "meters"
    v = empty.createVariable("y_reference_point", "f8")
    v.long_name = "y of the master reference point at southeast corner of field"
    v.units = "meters"

    # lat/lon positional variables
    v = empty.createVariable("lat_img_ne", "f8")
    v.long_name = "Latitude of northeast corner of image"
    v.units = "degrees_north"
    v = empty.createVariable("lat_img_nw", "f8")
    v.long_name = "Latitude of northwest corner of image"
    v.units = "degrees_north"
    v = empty.createVariable("lat_img_se", "f8")
    v.long_name = "Latitude of southeast corner of image"
    v.units = "degrees_north"
    v = empty.createVariable("lat_img_sw", "f8")
    v.long_name = "Latitude of southwest corner of image"
    v.units = "degrees_north"
    v = empty.createVariable("lat_reference_point", "f8")
    v.long_name = "Latitude of the master reference point at southeast corner of field"
    v.units = "degrees_north"
    v = empty.createVariable("lon_img_ne", "f8")
    v.long_name = "Longitude of northeast corner of image"
    v.units = "degrees_east"
    v = empty.createVariable("lon_img_nw", "f8")
    v.long_name = "Longitude of northwest corner of image"
    v.units = "degrees_east"
    v = empty.createVariable("lon_img_se", "f8")
    v.long_name = "Longitude of southeast corner of image"
    v.units = "degrees_east"
    v = empty.createVariable("lon_img_sw", "f8")
    v.long_name = "Longitude of southwest corner of image"
    v.units = "degrees_east"
    v = empty.createVariable("lon_reference_point", "f8")
    v.long_name = "Longitude of the master reference point at southeast corner of field"
    v.units = "degrees_east"

    # add global attributes
    empty.title = "None given"
    empty.created_by = "nobody"
    empty.history = "created with Python"

    empty.close()

# extract spectral profiles from environmentlogger.json
def irradiance_time_extractor(camera_type, envlog_file):
    # For the environmental logger records after 04/26/2016, there would be 24 files per day (1 file per hour, 5 seconds per record)
    # Convert json fiel to dictionary format file
    with open(envlog_file, "r") as fp:
        lines = fp.readlines()
        slines = "".join(lines)
        js = json.loads(slines)

    # assume that time stamp follows in 5 second increments across records since 5 sec/record
    num_readings = len(js["environment_sensor_readings"])
    if "spectrometers" in js["environment_sensor_readings"][0]:
        if camera_type == "swir_new":
            num_bands = len(js["environment_sensor_readings"][0]["spectrometers"]["NIRQuest-512"]["spectrum"])
        else:
            num_bands = len(js["environment_sensor_readings"][0]["spectrometers"]["FLAME-T"]["spectrum"])
    else:
        num_bands = len(js["environment_sensor_readings"][0]["spectrometer"]["spectrum"])

    spectra = np.zeros((num_readings, num_bands))
    times = []
    for idx in range(num_readings):
        # read time stamp
        time_current = js["environment_sensor_readings"][idx]["timestamp"]
        C = time_current.replace("."," ").replace("-"," ").replace(":","")
        ArrayTime=C.split(" ")
        time_current_r = int(ArrayTime[3])
        times.append(time_current_r)

        # read spectrum from irridiance sensors
        if "spectrometers" in js["environment_sensor_readings"][idx]:
            if camera_type == "swir_new":
                spectrum = js["environment_sensor_readings"][0]["spectrometers"]["NIRQuest-512"]["spectrum"]
            else:
                spectrum = js["environment_sensor_readings"][idx]["spectrometers"]["FLAME-T"]["spectrum"]
        else:
            spectrum = js["environment_sensor_readings"][idx]["spectrometer"]["spectrum"]

        spectra[idx,:] = spectrum

    return times, spectra

# replace rfl_img variable in netcdf with given matrix
def update_netcdf(inp, rfl_data, camera_type):
    print("Updating %s" % inp)

    out = inp.replace(".nc", "_newrfl.nc")

    with Dataset(inp) as src, Dataset(out, "w") as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            if name == "Google_Map_View":
                continue

            # Create variables
            var_dict = (src[name].__dict__)
            if '_FillValue' in var_dict.keys():
                x = dst.createVariable(name, variable.datatype, variable.dimensions, fill_value=var_dict['_FillValue'])
                del var_dict['_FillValue']
            else:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)

            # Set variables to values
            if name != "rfl_img":
                print("...%s" % name)
                dst[name][:] = src[name][:]
            else:
                if camera_type=='vnir_old':
                    print("...%s (subset)" % name)
                    dst[name][:679,:,:] = rfl_data
                    # 679-955 set to NaN
                    print("...NaNs")
                    dst[name][679:,:,:] = np.nan

                elif camera_type == "vnir_middle":
                    print("...%s (subset)" % name)
                    dst[name][:662,:,:] = rfl_data
                    # 679-955 set to NaN
                    print("...NaNs")
                    dst[name][662:,:,:] = np.nan
                else:
                    print("...%s" % name)
                    dst[name][:] = rfl_data

            # copy variable attributes all at once via dictionary
            dst[name].setncatts(var_dict)

        if 'rfl_img' not in src.variables:
            print("...adding rfl_img")
            dst.createVariable("rfl_img", "f4")
            dst.variables['rfl_img'] = rfl_data

# calculate frame times and solar zenith angle
def prepare_header_data(hdr_file, dataset_date):
    date_tuple = datetime.strptime(dataset_date, "%Y-%m-%d").timetuple()
    unix_basetime = date(year=1970, month=1, day=1)
    time_split = date(year=date_tuple.tm_year, month=date_tuple.tm_mon, day=date_tuple.tm_mday) - unix_basetime

    # Extract time of each frame from frameIndex
    framelist = []
    frame_file = hdr_file.replace("_raw.hdr", "_frameIndex.txt")
    with open(frame_file, 'r') as frames:
        for fl in frames.readlines()[1:]:
            hour_tuple = datetime.strptime(fl.split()[1], "%H:%M:%S").timetuple()
            framelist.append((time_split.total_seconds() + hour_tuple.tm_hour*3600.0 + hour_tuple.tm_min*60.0 +
                              hour_tuple.tm_sec) / (3600.0*24.0))

    sza = [solar_zenith_angle(datetime(year=1970, month=1, day=1) + timedelta(days=ftime)) for ftime in framelist]

    return {
        'frametime': framelist,
        'solar_zenith_angle': sza
    }

# populate empty netCDF headers with metadata & geo information
def update_netcdf_headers(nc_file, geodata, header_data):
    with Dataset(nc_file, 'a', mmap=False) as src:
        src["lat_img_ne"][...] = geodata["bbox_geojson"]["coordinates"][0][0][1]
        src["lon_img_ne"][...] = geodata["bbox_geojson"]["coordinates"][0][0][0]
        src["lat_img_nw"][...] = geodata["bbox_geojson"]["coordinates"][0][1][1]
        src["lon_img_nw"][...] = geodata["bbox_geojson"]["coordinates"][0][1][0]
        src["lat_img_sw"][...] = geodata["bbox_geojson"]["coordinates"][0][2][1]
        src["lon_img_sw"][...] = geodata["bbox_geojson"]["coordinates"][0][2][0]
        src["lat_img_se"][...] = geodata["bbox_geojson"]["coordinates"][0][3][1]
        src["lon_img_se"][...] = geodata["bbox_geojson"]["coordinates"][0][3][0]

        src["x_img_ne"][...] = geodata["x_coordinates"][-1]
        src["y_img_ne"][...] = geodata["y_coordinates"][0]
        src["x_img_nw"][...] = geodata["x_coordinates"][0]
        src["y_img_nw"][...] = geodata["y_coordinates"][0]
        src["x_img_sw"][...] = geodata["x_coordinates"][0]
        src["y_img_sw"][...] = geodata["y_coordinates"][-1]
        src["x_img_se"][...] = geodata["x_coordinates"][-1]
        src["y_img_se"][...] = geodata["y_coordinates"][-1]

        src["x_pxl_sz"][...] = geodata["x_pixel_size"]
        src["y_pxl_sz"][...] = geodata["y_pixel_size"]
        src["x_reference_point"][...] = geodata["x_reference"]
        src["y_reference_point"][...] = geodata["y_reference"]

        # TODO: hyperspectral_calculation has x/y swapped for these
        src["x"][...] = geodata["y_coordinates"]
        src["y"][...] = geodata["x_coordinates"]

        src["frametime"][...] = header_data['frametime']
        src["solar_zenith_angle"][...] = header_data['solar_zenith_angle']

# replace rfl_img variable in netcdf with given matrix
def update_netcdf_band(inp, band, band_data, camera_type):
    if band % 10 == 0:
        print("band: %s" % band)

    if camera_type=='vnir_old':
        inp["rfl_img"][band,:,:] = band_data
    elif camera_type == "vnir_middle":
        inp["rfl_img"][band,:,:] = band_data
    else:
        inp["rfl_img"][band,:,:] = band_data

# Get the rgb bands from netcdf file and save as jpg image for quick view
def convert_netcdf_to_jpg(out_file, out_path, camera_type, timestamp):
    out_rgb = os.path.join(out_path, "%s_netcdf_L1_ua-mac_%s.jpg" % (camera_type.split("_")[0], timestamp)) # define output path and image name

    # extract corresponding rgb bands from different netcdf files (VNIR or SWIR)
    f = Dataset(out_file)
    img_arr = np.asarray(f.variables['rfl_img'][:])
    rgb = np.zeros((img_arr.shape[1],img_arr.shape[2],3),dtype=np.uint8)

    if camera_type.startswith("swir"):
        r = img_arr[24,:,:]; r *= 255.0/r.max()
        g = img_arr[51,:,:]; g *= 255.0/g.max()
        b = img_arr[120,:,:]; b *= 255.0/b.max()
    else:
        r = img_arr[376,:,:]; r *= 255.0/r.max()
        g = img_arr[235,:,:]; g *= 255.0/g.max()
        b = img_arr[95,:,:]; b *= 255.0/b.max()

    rgb[:,:,0] = r; rgb[:,:,1] = g; rgb[:,:,2] = b

    out_img = Image.fromarray(rgb)
    out_img.save(out_rgb)

    # free up memory
    del img_arr; r; g; b; rgb

# apply calibration algorithm to the raw data
def apply_calibration(raw_filepath, out_file, metadata=None):
    if not os.path.exists(out_file):
        print("Error: %s does not exist to calibrate" % out_file)
        return

    print("Calibrating %s" % raw_filepath)

    # get necessary paths from path to _raw file
    raw_dir = os.path.dirname(raw_filepath)
    raw_file = os.path.basename(raw_filepath)
    md_file = os.path.join(raw_dir, "%s_metadata.json" % raw_file[:-4])
    date = raw_filepath.split("/")[-3]
    envlog_dir = os.path.join(raw_root, "EnvironmentLogger/%s" % date)

    # determine type of sensor and age of camera
    if raw_filepath.find("VNIR") > -1:
        if date < "2018-08-18":
            camera_type = "vnir_old"
            num_spectral_bands = 955
            num_bands_irradiance = 1024
            image_scanning_time   =  540
        elif "2018-08-18" <= date < "2019-02-26":
            camera_type = "vnir_middle"
            num_spectral_bands = 939
            num_bands_irradiance = 1024
            image_scanning_time   =  540
        else:
            camera_type = "vnir_new"
            num_spectral_bands = 939
            num_bands_irradiance = 3648
            # it is obverved that it takes an average of 3.5 mins/scan  = 210 seconds
            image_scanning_time   =  210
    else:
        if date < "2019-02-26": # Note that no calibration models are available for old&middle swir data
            camera_type = "swir_old_middle"
        else:
            camera_type = "swir_new"
            num_spectral_bands = 275
            num_bands_irradiance = 512
            image_scanning_time   =  210

    print("MODE: ---------- %s ----------" % camera_type)

    # load the raw data set
    print("Loading %s.hdr" % raw_filepath)
    hdr_path = raw_filepath +'.hdr'
    try:
        raw = envi.open(hdr_path)
        img_DN = raw.open_memmap()
    except IOError:
        print('No such file named %s' % raw_filepath)

    if metadata is None:
        # TEST CODE
        json_file = raw_filepath.replace("_raw", "_metadata.json")
        geo = test_pixel2Geographic(json_file, hdr_path, camera_type.split("_")[0].upper())
    else:
        geo = pixel2Geographic(metadata, hdr_path, camera_type.split("_")[0].upper())

    print("Updating netCDF headers")
    header_data = prepare_header_data(hdr_path, date)
    update_netcdf_headers(out_file, geo, header_data)

    # Since no calibration models are available for swir_old_middle, directly convert the raw data to netcdf
    if camera_type =="swir_old_middle":
        img_DN = np.rollaxis(img_DN, 2, 0)
        # Generate output path and call the netCDF conversion function, convert the raw old&middle swir data to netcdf
        update_netcdf(out_file, img_DN, camera_type)
        # free up memory
        del img_DN

    else: # apply pre-computed calibration models
        # Load the previously created calibration models based on the camera_type
        best_matched = os.path.join(calib_root + "/" + "calibration_new", camera_type, 'best_matched_index.npy')
        bias = os.path.join(calib_root + "/" + "calibration_new", camera_type, 'bias_coeff.npy')
        gain = os.path.join(calib_root + "/" + "calibration_new", camera_type, 'gain_coeff.npy')
        # read EnvLog data
        print("Reading EnvLog files in %s" % envlog_dir)
        envlog_tot_time =  []
        envlog_spectra = np.array([], dtype=np.int64).reshape(0, num_bands_irradiance)
        for ef in os.listdir(envlog_dir):
            if ef.endswith("environmentlogger.json"):
                time, spectrum = irradiance_time_extractor(camera_type,os.path.join(envlog_dir, ef))
                envlog_tot_time += time
                # print("concatenating %s onto %s" % (spectrum.shape, envlog_spectra.shape))
                envlog_spectra = np.vstack([envlog_spectra, spectrum])

        # Find the best match time range between image time stamp and EnvLog time stamp
        num_irridiance_record = int(image_scanning_time/5)   # 210/5=4.2  ---->  5 seconds per record

        # concatenation of hour mins and seconds of the image time stamp (eg., 12-38-49 to 123849)
        with open(md_file) as json_file:
            img_meta_data = json.load(json_file)
        meta_data_time = img_meta_data['lemnatec_measurement_metadata']['gantry_system_variable_metadata']['time']
        image_time     = meta_data_time[-8:]
        image_time     = int(image_time.replace(":",""))

        # compute the absolute difference between
        print("Computing mean spectrum")
        abs_diff_time = np.zeros((len(envlog_tot_time)))
        for k in range(len(envlog_tot_time)):
            abs_diff_time[k] = abs(image_time  - envlog_tot_time[k])
        ind_closet_time = np.argmin(abs_diff_time)  # closest time index
        mean_spectrum  = np.mean(envlog_spectra[ind_closet_time : ind_closet_time + num_irridiance_record-1, :], axis=0)

        # load pre-computed the best matched index between image and irradiance sensor spectral bands
        best_matched_index = np.load(best_matched)
        test_irridance     =  mean_spectrum[best_matched_index.astype(int).tolist()]
        test_irridance_re  = np.resize(test_irridance, (1, num_spectral_bands))

        # load and apply precomputed coefficient to convert irradiance to DN
        b = np.load(bias)
        g = np.load(gain)
        if camera_type == "vnir_old":
            test_irridance_re = test_irridance_re[:,0:679]
            img_DN = img_DN[:,:,0:679]
        if camera_type == "vnir_middle":
            test_irridance_re = test_irridance_re[:,0:662]
            img_DN = img_DN[:,:,0:662]

        irrad2DN = (g * test_irridance_re) + b

        # reflectance computation
        print("Computing reflectance and updating %s" % out_file)
        with Dataset(out_file, 'a', mmap=False) as src:
            # TODO: If file size is small enough, don't do this per-band but all at once
            for band_ind in range(img_DN.shape[2]):
                calibrated_band = img_DN[:,:,band_ind] / irrad2DN[:,band_ind]
                update_netcdf_band(src, band_ind, calibrated_band, camera_type)

        # Generate jpg format quick view image of final netcdf file
        # TODO: Runs out of memory on large files
        # convert_netcdf_to_jpg(out_file,out_path,camera_type,timestamp)

        del img_DN
        print("All done")
