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


# extract spectral profiles from environmentlogger.json
def irradiance_time_extractor(envlog_file):
    # For the environmental logger records after 04/26/2016, there would be 24 files per day (1 file per hour, 5 seconds per record)
    with open(envlog_file, "r") as fp:
        lines = fp.readlines()
        slines = "".join(lines)
        js = json.loads(slines)

    # assume that time stamp follows in 5 second increments across records since 5 sec/record
    num_readings = len(js["environment_sensor_readings"])
    if "spectrometers" in js["environment_sensor_readings"][0]:
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
                spectrum = js["environment_sensor_readings"][idx]["spectrometers"]["FLAME-T"]["spectrum"]
            else:
                spectrum = js["environment_sensor_readings"][idx]["spectrometer"]["spectrum"]
            spectra[idx,:] = spectrum
    
    return times, spectra

# apply calibration algorithm to the raw data
def apply_calibration(raw_filepath):
    print("Applying calibration to %s" % raw_filepath)
    # get necessary paths from path to _raw file
    raw_dir = os.path.dirname(input_path)
    raw_file = os.path.basename(input_path)
    md_file = os.path.join(raw_dir, "%s_metadata.json" % raw_file[:-4])
    date = raw_filepath.split("/")[-3]
    envlog_dir   = "/home/extractor/sites/ua-mac/raw_data/EnvironmentLogger/%s" % date

    # determine type of sensor and age of camera
    if raw_filepath.find("VNIR") > -1:
        if date < "2018-08-18":
            camera_type = "vnir_old"
            num_spectral_bands = 955
        else:
            camera_type = "vnir_new"
            num_spectral_bands = 939
    else:
        if date < "2018-08-17":
            camera_type = "swir_new" # TODO: does swir_old need separate calibration?
            num_spectral_bands = 273
        else:
            camera_type = "swir_new"
            num_spectral_bands = 275

    best_matched = os.path.join("calibration_new", camera_type, 'best_matched_index.npy')
    bias = os.path.join("calibration_new", camera_type, 'bias_coeff.npy')
    gain = os.path.join("calibration_new", camera_type, 'gain_coeff.npy')

    # load the raw data set
    print("Loading %s.hdr" % raw_filepath)
    try:
        raw = envi.open(raw_filepath +'.hdr')
        img_DN  = raw.load()
        #head_file = envi.read_envi_header(data_fullpath +'.hdr')
    except IOError:
        print('No such file named %s' % raw_filepath)

    # read EnvLog data
    print("Reading EnvLog files in %s" % envlog_dir)
    # num_bands_irradiance = 3648 # in new irradiance sensor, there are 3648 spectral bands
    envlog_tot_time =  []
    envlog_spectra = None
    for ef in os.listdir(envlog_dir):
        if ef.endswith("environmentlogger.json"):
            time, spectrum = irradiance_time_extractor(os.path.join(envlog_dir, ef))
            envlog_tot_time += time
            if envlog_spectra is None:
                envlog_spectra = np.array([], dtype=np.int64).reshape(0, spectrum.shape[1])
            envlog_spectra = np.vstack([envlog_spectra, spectrum])

    # Find the best match time range between image time stamp and EnvLog time stamp
    image_scanning_time   =  210 # it is obverved that it takes an average of 3.5 mins/scan  = 210 seconds
    num_irridiance_record = int(image_scanning_time/5)   # 210/5=4.2  ---->  5 seconds per record

    # concatenation of hour mins and seconds of the image time stamp (eg., 12-38-49 to 123849)
    with open(md_file) as json_file:
        img_meta_data = json.load(json_file)
    meta_data_time = img_meta_data['lemnatec_measurement_metadata']['gantry_system_variable_metadata']['time']
    image_time     = meta_data_time[-8:]
    image_time     = int(image_time.replace(":",""))

    # computer the absolute difference between
    abs_diff_time = np.zeros((len(envlog_tot_time)))
    for k in range(len(envlog_tot_time)):
        abs_diff_time[k] = abs(image_time  - envlog_tot_time[k])
    ind_closet_time = np.argmin(abs_diff_time)  # closest time index
    mean_spectrum   = np.mean(envlog_spectra[ind_closet_time : ind_closet_time + num_irridiance_record-1, :], axis=0)

    # load pre-computed the best matched index between image and irradiance sensor spectral bands
    best_matched_index = np.load(best_matched)
    test_irridance     =  mean_spectrum[best_matched_index.astype(int).tolist()]
    test_irridance_re  = np.resize(test_irridance, (1, num_spectral_bands))

    # load and apply precomputed coefficient to convert irradiance to DN
    b = np.load(bias)
    g = np.load(gain)
    irrad2DN = g*test_irridance_re + b

    # reflectance computation
    Ref  = (img_DN)/(irrad2DN) # Ref stands for the reflectance image

    # prepare output paths
    out_path = os.path.dirname(raw_filepath.replace("raw_data", "Level_1").replace("SWIR", "swir_netcdf").replace("VNIR", "vnir_netcdf"))
    out_file = os.path.join(out_path, 'ref_%s.npy' % raw_file) # os.path.join('ref_%s.hdr' % raw_file)
    if not os.path.isdir(out_path):
        os.makedirs(out_path)

    # TODO: Write to nc file
    """from ncdump -h *.nc

    variables:
        float rfl_img(wavelength, y, x) ;
            rfl_img:long_name = "Reflectance of image" ;
            rfl_img:standard_name = "surface_albedo" ;
            rfl_img:units = "1" ;


    with Dataset(rfl_image, "w") as netCDF_handle:
        netCDF_handle.write(Ref)
    """

    # save as ENVI file (RGB bands: 392, 252, 127)
    # envi.save_image(out_file, Ref, dtype=np.float32, interleave='bil', force = 'True', metadata=head_file)

    # Save Ref as a .npy file
    np.save(out_file, Ref)


# TODO: This will come from the extractor message
input_path  = "/home/extractor/sites/ua-mac/raw_data/SWIR/2017-04-16/2017-04-16__11-50-46-707/c6079666-b686-4481-9a4f-0663f5f43a6a_raw"
# "/home/extractor/sites/ua-mac/raw_data/SWIR/2018-08-17/2018-08-17__15-44-34-088/c283d6ef-9ccf-4f8f-b291-bb753e50c0ac_raw"
# "/home/extractor/sites/ua-mac/raw_data/VNIR/2018-01-23/2018-01-23__06-59-08-955/757ceefe-afc6-4b80-b126-e7b7abf85a4d_raw"
# "/home/extractor/sites/ua-mac/raw_data/VNIR/2018-08-18/2018-08-18__11-11-41-890/c5f4d50f-44ad-4e23-9d92-f10e62110ac7_raw"
apply_calibration(input_path)
