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


raw_root = "/home/extractor/hs_calib"
calib_root = "/home/extractor"


def create_empty_netCDF(raw_file):
    out_path = os.path.dirname(raw_file.replace("raw_data", "Level_1").replace("SWIR", "swir_netcdf").replace("VNIR", "vnir_netcdf"))
    timestamp = raw_file.split("/")[-2]
    if raw_file.find("VNIR") > -1:
        camera_type = "vnir"
    else:
        camera_type = "swir"
    out_file = os.path.join(out_path, "%s_netcdf_L1_ua-mac_%s.nc" % (camera_type, timestamp))

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

    # add dimensions
    empty.createDimension("wavelength", hdr_bands)
    empty.createDimension("x", hdr_samples)
    empty.createDimension("y", hdr_lines)

    # add Variable double temperature(lat,lon)
    t = empty.createVariable("xps_img", "u2", ("wavelength", "x", "y"))
    t.long_name = "Exposure counts"
    t.meaning = "Counts on scale from 0 to 2^16-1 = 65535"
    t.units = "1"

    empty.createVariable("rfl_img", "f4", ("wavelength", "x", "y"))

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


# apply calibration algorithm to the raw data
def apply_calibration(raw_filepath):
    print("Calibrating %s" % raw_filepath)

    # get necessary paths from path to _raw file
    raw_dir = os.path.dirname(raw_filepath)
    raw_file = os.path.basename(raw_filepath)
    md_file = os.path.join(raw_dir, "%s_metadata.json" % raw_file[:-4])
    date = raw_filepath.split("/")[-3]
    timestamp = raw_filepath.split("/")[-2]
    envlog_dir   = os.path.join(raw_root, "EnvironmentLogger/%s" % date)

    #     prepare output paths
    print("Generating output")
    out_path = os.path.dirname(raw_filepath.replace("raw_data", "Level_1").replace("SWIR", "swir_netcdf").replace("VNIR", "vnir_netcdf"))
    if not os.path.isdir(out_path):
        os.makedirs(out_path)
    if os.path.isfile(out_path):
        print("Output file already exists: skipping "+out_path)
        return

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
    try:
        raw = envi.open(raw_filepath +'.hdr')
        #        img_DN  = raw.load()
        img_DN = raw.open_memmap()
        #head_file = envi.read_envi_header(data_fullpath +'.hdr')
    except IOError:
        print('No such file named %s' % raw_filepath)

    # Apply calibration procedure if camera_type == vnir_old, vnir_middle, vnir_new or swir_new. Since no calibration models are available for
    # swir_old and swir_middle, so directly convert old&middel SWIR raw data to netcdf format
    if camera_type =="swir_old_middle":
        # Convert the raw swir_old and swir_middle data to netCDF
        img_DN = np.rollaxis(img_DN, 2, 0)
        # Generate output path and call the netCDF conversion function, convert the raw old&middle swir data to netcdf
        out_path = os.path.dirname(raw_filepath.replace("raw_data", "Level_1").replace("SWIR", "swir_netcdf").replace("VNIR", "vnir_netcdf"))
        out_file = os.path.join(out_path, "%s_netcdf_L1_ua-mac_%s.nc" % (camera_type.split("_")[0], timestamp))
        update_netcdf(out_file, img_DN, camera_type)

        # free up memory
        del img_DN

    else: # when camera_type == vnir_old, vnir_middle, vnir_new or swir_new, apply pre-computed calibration models
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
        out_file = os.path.join(out_path, "%s_netcdf_L1_ua-mac_%s.nc" % (camera_type.split("_")[0], timestamp))
        print("Computing reflectance and updating %s" % out_file)
        with Dataset(out_file, 'a', mmap=False) as src:
            for band_ind in range(img_DN.shape[2]):
                calibrated_band = img_DN[:,:,band_ind] / irrad2DN[:,band_ind]
                calibrated_band = np.swapaxes(calibrated_band, 0, 1)
                update_netcdf_band(src, band_ind, calibrated_band, camera_type)

        print("All done")

        # save as ENVI file (RGB bands: 392, 252, 127)
        # VNIR - 376 (R) 235 (G) 95 (B)
        out_file = os.path.join('ref_%s.hdr' % raw_file)
        envi.save_image(out_file, Ref, dtype=np.float32, interleave='bil', force = 'True', metadata=head_file)


if __name__ == "__main__":
    # TODO: This will come from the extractor message
    input_paths = [
        os.path.join(raw_root, "VNIR/2019-05-23/2019-05-23__13-52-42-629/47265f0b-7a51-43a3-8f7a-bdf11bd6fe38_raw"),
        #os.path.join(raw_root, "VNIR/2018-08-18/2018-08-18__11-11-41-890/c5f4d50f-44ad-4e23-9d92-f10e62110ac7_raw"),
    ]

    for p in input_paths:
        create_empty_netCDF(p)
        apply_calibration(p)

