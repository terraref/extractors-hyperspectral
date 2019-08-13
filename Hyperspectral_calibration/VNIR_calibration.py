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
from numpy import loadtxt
from scipy.interpolate import interp1d

#%% Notation
'''
data_folder                : raw image data directory
fileName                   : the raw image name in the data folder
data_fullpath              : image directory + image name
Environmental_Loggar_files : a list of available Environmental Logga files
Environmental_Loggar_dir   : the directory of environmental Loggar files that store irradance values
irradiance_time_extractor  : a function to extract spectral profiles from environmentlogger.json

'''
#%% set the raw data directory
data_folder     = r"I:\terra_ref\Datasets\caibration_data_05_2019\2019-02-28__14-32-24-455"
fileName        = '57b5be99-e866-4f05-85c9-f35d66dc392b_raw'
data_fullpath   = data_folder + '\\' + fileName 

#  Environmental_logger file for raw VNIR data
Environmental_Loggar_files = os.listdir(r"I:\terra_ref\Datasets\caibration_data_05_2019\2019-02-28_Environmental_Loggar")
Environmental_Loggar_dir   = r"I:\terra_ref\Datasets\caibration_data_05_2019\2019-02-28_Environmental_Loggar/"
save_name                  = 'ref_' + fileName + '.hdr'
# TODO: should this only be run on new VNIR sensor given these bands?
num_spectral_bands_VNIR    = 939  # in new VNIR sensor, there are 939 spectral bands
num_bands_irradiance       = 3648 # in new irradiance sensor, there are 3648 spectral bands

#%%
#For the environmental logger records after 04/26/2016, there would be 24 files per day (1 file per hour, 5 seconds per record)
# EnvironmentLogger info extraction

def  irradiance_time_extractor(FileName):
    fp = open(FileName, "r")
    lines = fp.readlines()
    slines = "".join(lines)
    js = json.loads(slines)
    #fp.close()
    
    firstTime=str()
    firstTime=js["environment_sensor_readings"][0]["timestamp"]
    firstTime=firstTime.replace("."," ").replace("-"," ")
    
    sz = len(js["environment_sensor_readings"])
    
    # assume that time stamp follows in 5 second increments across records since 5 sec/record
    num_bands = len(js["environment_sensor_readings"][0]["spectrometers"]["FLAME-T"]["spectrum"])
    nSpectrum = np.zeros((sz, num_bands))
    nTime     = []
    for idx in range(sz):
            # read time stamp
            time_current = js["environment_sensor_readings"][idx]["timestamp"]
            C = time_current.replace("."," ").replace("-"," ").replace(":","")
            ArrayTime=C.split(" ")
            time_current_r = int(ArrayTime[3])
            nTime.append(time_current_r)
            
            # read spectrum from irridiance sensors
            spectrum = js["environment_sensor_readings"][idx]["spectrometers"]["FLAME-T"]["spectrum"]
            nSpectrum[idx,:] = spectrum
    
    return nTime, nSpectrum
   
#%%
#  ***Apply calibration algorithm to the raw VNIR data ***

import spectral.io.envi as envi

# load the raw data set
try:
    raw = envi.open(data_fullpath +'.hdr')
    img_DN  = raw.load()
    head_file = envi.read_envi_header(data_fullpath +'.hdr')
    print('Load data successfully ...')
except IOError:
    print('No such file named %s', fileName)

#%% read enviromental loggar data
ntotal_time_Enviormental_logger  =  []
ntotal_Spectrum_Enviormental_logger  =  np.array([], dtype=np.int64).reshape(0, num_bands_irradiance)
for j in range(len(Environmental_Loggar_files)):
    FileName         = Environmental_Loggar_dir + Environmental_Loggar_files[j]
    nTime, nSpectrum = irradiance_time_extractor(FileName)
    ntotal_time_Enviormental_logger     = ntotal_time_Enviormental_logger + nTime
    ntotal_Spectrum_Enviormental_logger = np.vstack([ntotal_Spectrum_Enviormental_logger, nSpectrum])    
    
# Find the best match time range between image time stamp and Enviormental_logger time stamp
image_scanning_time   =  210 # it is obverved that it takes an average of 3.5 mins/scan  = 210 seconds  
num_irridiance_record = int(image_scanning_time/5)   # 210/5=4.2  ---->  5 seconds per record
#%
# concatenation of  hour  mins and seconds of the image time stamp (eg., 12-38-49 to 123849)
meta_data_name = data_folder + '\\' + fileName[:-4] + "_metadata.json"
with open(meta_data_name) as json_file:  
    img_meta_data = json.load(json_file)
meta_data_time = img_meta_data['lemnatec_measurement_metadata']['gantry_system_variable_metadata']['time']
image_time     = meta_data_time[-8:]
image_time     = int(image_time.replace(":",""))

# computer the absolute difference between 
abs_diff_time = np.zeros((len(ntotal_time_Enviormental_logger)))
for k in range(len(ntotal_time_Enviormental_logger)):
    abs_diff_time[k] = abs(image_time  - ntotal_time_Enviormental_logger[k])

ind_closet_time = np.argmin(abs_diff_time)  # closest time index
mean_spectrum   = np.mean(ntotal_Spectrum_Enviormental_logger[ind_closet_time : ind_closet_time + num_irridiance_record-1, :], axis=0)
#
best_matched_index = np.load('best_matched_index.npy') # load pre-computed the best matched index between image and irradiance sensor spectral bands
test_irridance     =  mean_spectrum[best_matched_index.astype(int).tolist()]
test_irridance_re  = np.resize(test_irridance, (1, num_spectral_bands_VNIR))

# load and apply precomputed coefficient to convert irradiance to DN
b = np.load('bias_coeff.npy') # load the bias coefficient
g = np.load('gain_coeff.npy') # load the gain coefficient
irrad2DN = g*test_irridance_re + b

# reflectance computation
Ref  = (img_DN)/(irrad2DN)#   Ref stands for the reflectance image
#
# save as ENVI file (RGB bands: 392, 252, 127) 
envi.save_image(data_folder + "\\" + save_name, Ref, dtype=np.float32, interleave='bil', force = 'True', metadata=head_file)

    
    
    