import json
import os
import numpy as np
import spectral.io.envi as envi
from PIL import Image
from netCDF4 import Dataset
# crop to plot helper lib
from terra_common import CoordinateConverter as CC
# libs for ML hs reduction
from sklearn.utils.extmath import randomized_svd
import math
import h5py

# this was developed as a team effort between:
#   Michele Cosi - SLU (Vasit Sagan group)
#   Maxwell Burnette - NCSA
#   Abby Stylianou - SLU
#   Robert Pless - GWU
#   Roxana Leontie - GWU

# paths:
raw_root = "/home/extractor/hs_calib"
calib_root = "/home/extractor/"
out_root = "/home/extractor/hs_calib/S9_ProcessTest"
env_log = "/home/extractor/hs_calib"
# variable to determine the dimensionality reduction of the output calibrated data
nc_reduction = 10
# terra common coord converter used by find_crop_position
cc = CC()

def create_empty_ncdf(out_file, x, y, bands):
    empty = Dataset(out_file, "w", format="NETCDF4")
    # add dimensions
    empty.createDimension("wavelength", bands)
    empty.createDimension("x", math.ceil(x/nc_reduction))
    empty.createDimension("y", math.ceil(y/nc_reduction))

    empty.createVariable("rfl_img", "f4", ("wavelength", "x", "y"))

    # add global attributes
    empty.title = "VNIR plot level calibrated file"
    empty.created_by = "NCSA & SLU & GWU"
    empty.history = "created with Python"
    empty.close()

def update_netcdf_band(out_file, band, band_data):
    # if band % 10 == 0:
    #     print("band: %s" % band)
    #     print('size band_data %s' % band_data.size)
    with Dataset(out_file, 'a', mmap=False) as src:
        src["rfl_img"][band, :, :] = band_data[::nc_reduction, ::nc_reduction]

def irradiance_time_extractor(camera_type, envlog_file):

    # extract spectral profiles from environmentlogger.json

    # For the environmental logger records after 04/26/2016, there would be 24 files per day (1 file per hour, 5 seconds per record)
    # Convert json field to dictionary format file
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
        C = time_current.replace(".", " ").replace("-", " ").replace(":", "")
        ArrayTime = C.split(" ")
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

        spectra[idx, :] = spectrum

    return times, spectra

def pixel_to_fieldCoord(raw_filepath):

    print("Computing field coordinates for %s" % raw_filepath)

    # get necessary paths from path to _raw file
    raw_dir = os.path.dirname(raw_filepath)
    raw_file = os.path.basename(raw_filepath)
    json_path = os.path.join(raw_dir, "%s_metadata.json" % raw_file[:-4])
    hdr_path = raw_filepath + '.hdr'
    timestamp = raw_filepath.split("/")[-2]

    # ----- GET VARIABLE METADATA -----
    with open(json_path) as json_meta:
        meta = json.loads(json_meta.read())["lemnatec_measurement_metadata"]

    x_gantry_pos, y_gantry_pos = None, None

    x_camera_pos = float(meta['sensor_fixed_metadata']["location in camera box x [m]"])
    y_camera_pos = float(meta['sensor_fixed_metadata']["location in camera box y [m]"])
    x_gantry_pos = float(meta['gantry_system_variable_metadata']['position x [m]'])
    y_gantry_pos = float(meta['gantry_system_variable_metadata']['position y [m]'])
    scan_dir = str(meta['gantry_system_variable_metadata']["scanDirectionIsPositive"])

    # ----- CALCULATE POSITIONS -----
    # equivalent to https://github.com/GWUvision/OPEN/blob/8e74c9ff438afc9a0fad8a73157b0fc5c53d1c3b/crop/stereo.py#L152 in stereo
    # this is the start position
    x_absolute_pos = x_gantry_pos + x_camera_pos
    y_absolute_pos = y_gantry_pos + y_camera_pos

    # print("x_camera_pos: ", x_camera_pos)
    # print("x_gantry_pos: ", x_gantry_pos)
    # print("y_camera_pos: ", y_camera_pos)
    # print("y_gantry_pos: ", y_gantry_pos)
    # print("x_absolute_pos center of the image: ", x_absolute_pos)
    # print("y_absolute_pos center of the image: ", y_absolute_pos)

    # ----- GET HDR METADATA -----
    with open(hdr_path) as hdr_data:
        overall = hdr_data.readlines()
    for members in overall:
        if "samples" in members:
            x_pixel_num = int(members.split("=")[-1].strip("\n"))
            #print("x_pixel_num = ", x_pixel_num)
        if "lines" in members:
            y_pixel_num = int(members.split("=")[-1].strip("\n"))
            #print("y_pixel_num = ", y_pixel_num)

    # ----- PREPARE METADATA -----
    if raw_filepath.find("SWIR") > -1:
        x_pixel_size = 1.930615052e-3
    else:  # VNIR
        # seasons 9: FOV is 1m at distance 2m (from the canopy) we have 1600 pixels for x 1/1600 = 0.000625
        # if the ground can be seen the FOV should be calculated as:
        # updated_fov = z_absolute_pos * 1 (meter) / 2 (meters)
        # x_pixel_size = updated_fov / x_pixel_num

        # since the camera is focused at 2m (canopy) the decision was made to assume a fov of 1 meter,
        # otherwise we would need to know when we have canopy closure in the season and use the 1 meter one only then
        x_pixel_size = 1 / x_pixel_num

    y_pixel_size = 0.98526434004512529576754637665e-3

    # Roxana's notes: for the first half of the image we substract the distance in pixels and for the last half we add them
    x_final_result = []
    x_half_px = x_pixel_num / 2
    for x in range(x_pixel_num):
        if x <= x_half_px:
            x_final_result_value = x_absolute_pos + ((-x + x_half_px) * -x_pixel_size)
        else:
            x_final_result_value = x_absolute_pos + ((x - x_half_px) * x_pixel_size)
        x_final_result.append(x_final_result_value)
    x_final_result = np.array(x_final_result)

    if scan_dir == "True":
        #print("scan direction: e to w")
        # calculate the position of each pixel along the scan direction in the field coordinates
        y_final_result = np.array([y * y_pixel_size for y in range(y_pixel_num)]) + y_absolute_pos + 0.16 # works for season 9 & 10
        # Roxana's notes: bounding box of the image in field coordinates for scanning e to w
        # first and last item in the list of coordinates for each measurement
        #print("x_final_result[0] ", x_final_result[0])
        #print("x_final_result[-1] ", x_final_result[-1])
        #print("y_final_result[0] ", y_final_result[0])
        #print("y_final_result[-1] ", y_final_result[-1])
    else:
        #print("scan direction: w to e")
        y_final_result = np.array([-y * y_pixel_size for y in range(y_pixel_num)]) + y_absolute_pos - 0.04 # works for season 9 & 10

        # Roxana's notes: bounding box of the image in field coordinates for scanning w to e
        # first and last item in the list of coordinates for each measurement
        #print("x_final_result[0] ", x_final_result[0])
        #print("x_final_result[-1] ", x_final_result[-1])
        #print("y_final_result[0] ", y_final_result[0])
        #print("y_final_result[-1] ", y_final_result[-1])

        #print("dim final x", x_final_result.shape)
        #print("dim final y", y_final_result.shape)

    return x_final_result, y_final_result

def max_range(seq, value):
    position = np.argwhere(seq == value)
    return position.min(), position.max() + 1

def generate_geojson(SE, SW, NE, NW):
    # Create a geojson file for viewing e.g. in QGIS
    return {
        "type": "Polygon",
        "coordinates": [[
            [NE[1], NE[0]],
            [NW[1], NW[0]],
            [SW[1], SW[0]],
            [SE[1], SE[0]],
            [NE[1], NE[0]]  # repeat final coordinate
        ]]
    }

def bounding_box_to_geographic(x_min, x_max, y_min, y_max, out_path):

    LATITUDE_TO_METER = 1 / (30.80716 * 3600)
    LONGITUDE_TO_METER = 1 / (25.906 * 3600)  # varies, but has been corrected based on the actual location of the field

    # hardcoded from the reference point shape file
    origin_x = 33.074543
    origin_y = -111.97479

    SE = (x_min * LATITUDE_TO_METER + origin_x, -(y_min) * LONGITUDE_TO_METER + origin_y)
    SW = (x_min * LATITUDE_TO_METER + origin_x, -(y_max) * LONGITUDE_TO_METER + origin_y)
    NE = (x_max * LATITUDE_TO_METER + origin_x, -(y_min) * LONGITUDE_TO_METER + origin_y)
    NW = (x_max * LATITUDE_TO_METER + origin_x, -(y_max) * LONGITUDE_TO_METER + origin_y)

    j_bbox = generate_geojson(SE, SW, NE, NW)

    with open(out_path, 'w') as shpfile:
        json.dump(j_bbox, shpfile)

def rotate_band(band, scan_dir):

    # raw data is stored aligned length x width we want it aligned width x length
    band = np.swapaxes(band, 0, 1)

    # we have to flip data vertically and horizontally for e -> w, and only vertically w -> e
    if scan_dir == "True":
        # print("scan dir e -> w")
        # flip vertically
        band = np.flipud(band)
        # flip horizontally
        band = np.fliplr(band)
    else:
        #print("scan dir w -> e")
        # flip vertically
        band = np.flipud(band)

    return band

def save_image(rgb, out_path, timestamp):

    rgb = np.asarray(rgb)
    # somehow the image creation swaps the axis. We want the image and data to be recoded with the same dimensions
    rgb = np.swapaxes(rgb, 0, 1)
    out_img = Image.fromarray(rgb)
    # saving the image as png takes 25x more space no need to do it because imgs are for visualization
    out_rgb = os.path.join(out_path, "%s.jpg" % timestamp)
    print("saving image: ", out_rgb)
    out_img.save(out_rgb)

def find_crop_position(raw_filepath):

    print("Finding plot crop positions")
    x_map, y_map = pixel_to_fieldCoord(raw_filepath)
    row_plot_num = np.zeros([y_map.size])
    plot_cols_range = []

    left_most_field_range, _ = cc.fieldPosition_to_fieldPartition(x_map[0], y_map[0])
    right_most_field_range, _ = cc.fieldPosition_to_fieldPartition(x_map[-1], y_map[-1])

    if left_most_field_range == right_most_field_range:
        plot_ranges = [left_most_field_range]
        plot_ranges_range = [[0, x_map.size - 1]]
    else:
        print("more than one range in the scan")
        col_plot_num = np.zeros([x_map.size])
        # Compute all ranges in the scan
        for i in range(x_map.size):
            plot_range, plot_col = cc.fieldPosition_to_fieldPartition(x_map[i], y_map[0])

        plot_ranges = np.unique(col_plot_num)
        plot_ranges = plot_ranges[plot_ranges != 0]
        plot_ranges_range = []
        for plot_range in plot_ranges:
            plot_ranges_range.append(max_range(col_plot_num, plot_range))

    # Compute all columns in the scan
    for i in range(y_map.size):
        plot_range, plot_col = cc.fieldPosition_to_fieldPartition(x_map[0], y_map[i])
        row_plot_num[i] = plot_col

    plot_cols = np.unique(row_plot_num)
    plot_cols = plot_cols[plot_cols != 0]
    for col in plot_cols:
        plot_cols_range.append(max_range(row_plot_num, col))
    crop_positions = {}

    # combine
    for plot_range, range_range in zip(plot_ranges, plot_ranges_range):
        for plot_col, col_range in zip(plot_cols, plot_cols_range):
            crop_positions[(plot_range, plot_col)] = [range_range, col_range]

    return crop_positions, x_map, y_map

def process_VNIR(raw_filepath, nc_calib, out_root_dir, raw_env_root):
    # this should work for TERRA Season 9 VNIR data
    # For a raw scan file given as a path input raw_filepath the output is:
    #   1. PCA reduction of the calibrate spectra at the plot level data as a h5
    #   2. Plot boundaries in field coordinates (json info),
    #   3. plot image RGB crop, and
    #   4. geojson for each plot?
    #   5. optional plot level calibrated nc file when nc_calib is true

    # ----- Read raw data -----
        # get necessary paths from path to _raw file

    raw_dir = os.path.dirname(raw_filepath)
    raw_file = os.path.basename(raw_filepath)
    md_file = os.path.join(raw_dir, "%s_metadata.json" % raw_file[:-4])
    date = raw_filepath.split("/")[-3]
    timestamp = raw_filepath.split("/")[-2]
    envlog_dir = os.path.join(raw_env_root, "EnvironmentLogger/%s" % date)


    # This is for season 9 VNIR
        # season 6 doesn't have any data -- cameras were not functional
    if raw_filepath.find("VNIR") > -1: # TODO: add date range here to make sure is for season 9?
        camera_type = "vnir_new"
        num_spectral_bands = 939
        num_bands_irradiance = 3648
        # it is obverved that it takes an average of 3.5 mins/scan  = 210 seconds
        # TODO: do we want to change this to be the exact number? data exist in metadata to compute?
        image_scanning_time = 210
    else:
        print("This file doesn't meet the requirements: VNIR season 9 date range: ", raw_filepath, date)
        return


    # print("MODE: ---------- %s ----------" % camera_type)

    # load the raw data set
    print("Loading %s.hdr" % raw_filepath)
    try:
        raw = envi.open(raw_filepath + '.hdr')
        # Roxana : this is where the image is loaded
        img_DN = raw.open_memmap()
    except IOError:
        print('No such file named %s' % raw_filepath)


    print(md_file)
    with open(md_file) as json_file:
        img_meta_data = json.load(json_file)

    # concatenation of hour mins and seconds of the image time stamp (eg., 12-38-49 to 123849)
    meta_data_time = img_meta_data['lemnatec_measurement_metadata']['gantry_system_variable_metadata']['time']
    # print(meta_data_time)
    image_time = meta_data_time[-8:]
    # print(image_time)
    image_time = int(image_time.replace(":", ""))
    # print(image_time)

    # keep track of scan direction e->w true, w->e false
    # this is useful for image rotations and offsets
    scan_dir = str(
        img_meta_data['lemnatec_measurement_metadata']['gantry_system_variable_metadata']['scanDirectionIsPositive'])

# ----- Determine the calibration parameters specific to this file -----

    # Load the previously created calibration models based on the camera_type
    best_matched = os.path.join(calib_root + "/" + "calibration_new", camera_type, 'best_matched_index.npy')
    bias = os.path.join(calib_root + "/" + "calibration_new", camera_type, 'bias_coeff.npy')
    gain = os.path.join(calib_root + "/" + "calibration_new", camera_type, 'gain_coeff.npy')

    # read EnvironmentLogger data
    print("Reading EnvLog files in %s" % envlog_dir)
    envlog_tot_time = []
    envlog_spectra = np.array([], dtype=np.int64).reshape(0, num_bands_irradiance)

    # TODO: this takes a long time with the purpose of computing the mean spectrum... maybe could be made faster?
    for ef in os.listdir(envlog_dir):
        if ef.endswith("environmentlogger.json"):
            time, spectrum = irradiance_time_extractor(camera_type, os.path.join(envlog_dir, ef))
            envlog_tot_time += time
            # print("concatenating %s onto %s" % (spectrum.shape, envlog_spectra.shape))
            envlog_spectra = np.vstack([envlog_spectra, spectrum])

    # Find the best match time range between image time stamp and EnvLog time stamp
    num_irridiance_record = int(image_scanning_time / 5)  # 210/5=4.2  ---->  5 seconds per record


    # compute the absolute difference between
    print("Computing mean spectrum")
    abs_diff_time = np.zeros((len(envlog_tot_time)))
    for k in range(len(envlog_tot_time)):
        abs_diff_time[k] = abs(image_time - envlog_tot_time[k])
    ind_closet_time = np.argmin(abs_diff_time)  # closest time index
    mean_spectrum = np.mean(envlog_spectra[ind_closet_time: ind_closet_time + num_irridiance_record - 1, :], axis=0)

    # load pre-computed the best matched index between image and irradiance sensor spectral bands
    best_matched_index = np.load(best_matched)
    test_irridance = mean_spectrum[best_matched_index.astype(int).tolist()]
    test_irridance_re = np.resize(test_irridance, (1, num_spectral_bands))

    # load and apply precomputed coefficient to convert irradiance to DN
    b = np.load(bias)
    g = np.load(gain)
    irrad2DN = (g * test_irridance_re) + b

# ----- Find the plot crop positions  -----
    crop_positions = {}
    crop_positions, x_map, y_map = find_crop_position(raw_filepath)

# ----- Hyperspectral PCA at the plot level -----

    for (plot_row, plot_col), (row_range, col_range) in crop_positions.items():
        # determine the plot number
        # plot_num = cc.fieldPartition_to_plotNum(plot_row, plot_col)
        # plot level the image array initialization
        rgb = np.zeros((img_DN.shape[1], col_range[1] - col_range[0], 3), dtype=np.uint8)

        # ----- Hyperspectral compression parameters -----

        # Lets find a multiplier to get down to about 40 mega-data-points.
        data_size = (col_range[1] - col_range[0], img_DN.shape[1], img_DN.shape[2])
        M = np.prod(data_size) / 40000000
        skipSize = math.ceil(math.sqrt(M))
        A = np.zeros((data_size[2], math.ceil(data_size[1] / skipSize), math.ceil(data_size[0] / skipSize)))

        for band_ind in range(img_DN.shape[2]):
            # apply calibration on the reduced size
            calibrated_band = img_DN[col_range[0]:col_range[1]:skipSize, ::skipSize, band_ind] / irrad2DN[:, band_ind]
            # reorient the data correctly based on the scan direction
            A[band_ind, :, :] = rotate_band(calibrated_band, scan_dir)

        At = np.transpose(A, [2, 1, 0])
        nC = A.shape[0]
        x = A.shape[1]
        y = A.shape[2]
        B = np.reshape(At, (x * y, nC), order="F")

        # the val returned here different than Matlab: V is already transposed and S is a row vector
        U, Sigma, Vt = randomized_svd(np.asmatrix(B), n_components=15, n_iter=5, random_state=None)
        # S needs a bit of work for multiplication... the result is an array
        S = np.diag(Sigma)
        reconstructedB = U @ S @ Vt  # equivalent to np.dot(np.dot(U,np.diag(S[::-1])),Vt)
        reconstructionError = B - reconstructedB
        maxError = reconstructionError.flatten().max()
        MSE = np.square(reconstructionError).mean(axis=None)
        VinvS = np.transpose(Vt) @ np.linalg.inv(S)

        # prepare output paths
        out_path = os.path.join(out_root_dir, "%s_%s" % (int(plot_row), int(plot_col)))

        if not os.path.isdir(out_path):
            os.makedirs(out_path)
        # calibrated plot level ncdf?
        if nc_calib:
            out_nc = os.path.join(out_path, "%s.nc" % timestamp)
            create_empty_ncdf(out_nc, img_DN.shape[1], col_range[1] - col_range[0], img_DN.shape[2])

        # apply calibration and PCA at the plot level
        UBig = 0
        for band_ind in range(img_DN.shape[2]):
            # ---- calibrate each band ----
            # use cols for rows and rows for cols because raw data is in reverse
            #calibrated_band = img_DN[col_range[0]:col_range[1], row_range[0]:row_range[1], band_ind] / irrad2DN[:, band_ind]
            calibrated_band = img_DN[col_range[0]:col_range[1], :, band_ind] / irrad2DN[:, band_ind]
            # reorient the data correctly based on the scan direction
            calibrated_band = rotate_band(calibrated_band, scan_dir)
            # save the calibrated band to the nc file
            if nc_calib:
                update_netcdf_band(out_nc, band_ind, calibrated_band)

            #print("calibrated band size ", calibrated_band.shape)
            # ---- dimensionality reduction ---
            Bslice_bigT = np.transpose(calibrated_band)
            VinvSrow = VinvS[band_ind, :]
            #print("VinvSrow size ", VinvSrow.shape)
            # Bslice_big.flatten('F') equivalent to Matlab Bslice_big(:)
            Utmp = np.transpose(np.asmatrix(Bslice_bigT.flatten('F'))) @ np.asmatrix(VinvSrow)
            UBig = UBig + Utmp
            # ---- save the RGB bands when we process them for image ----
            # is band Red?
            if band_ind == 376:
                # normalized
                rgb[:, :, 0] = calibrated_band * 255.0 / calibrated_band.max()
            # is band Green?
            if band_ind == 235:
                # normalized
                rgb[:, :, 1] = calibrated_band * 255.0 / calibrated_band.max()
            # is band Blue?
            if band_ind == 95:
                # normalized
                rgb[:, :, 2] = calibrated_band * 255.0 / calibrated_band.max()



        # ---- processing complete save the relevant output for the current plot ----
        save_image(rgb, out_path, timestamp)

        # should save UBig S V dataSize maxError MSE
        out_h5 = os.path.join(out_path, "%s.h5" % timestamp)
        h5f = h5py.File(out_h5, 'w')
        h5f.create_dataset('UBig', data=UBig, dtype='<f4')
        h5f.create_dataset('S', data=Sigma, dtype='<f4')
        h5f.create_dataset('Vt', data=Vt, dtype='<f4')
        h5f.create_dataset('maxError', data=maxError, dtype='<f4')
        h5f.create_dataset('MSE', data=MSE, dtype='<f4')
        h5f.create_dataset('dataSize', data=Bslice_bigT.shape)
        h5f.close()

        # write new json file

        # field coordinates (range) x ones should not change
        x_min = float(x_map[0])
        x_max = float(x_map[-1])
        y_box_1 = float(y_map[col_range[0]])
        y_box_2 = float(y_map[col_range[1]-1])


        if y_box_1 > y_box_2:
            y_min = y_box_2
            y_max = y_box_1
        else:
            y_min = y_box_1
            y_max = y_box_2

        add_metadata = {'field_box_meters': {'xmin': x_min, 'xmax': x_max, 'ymin': y_min, 'ymax': y_max}}
        # copy all the original json metadata in the plot json along with the plot boundary in field coordinates
        dst_json_data = dict(img_meta_data)
        dst_json_data.update(add_metadata)

        dst_json_path = os.path.join(out_path, "%s.json" % timestamp)
        with open(dst_json_path, 'w') as outfile:
            json.dump(dst_json_data, outfile, sort_keys=True, indent=4)

        # generate geojson
        out_gjson_path = os.path.join(out_path, "%s.geojson" % timestamp)

        bounding_box_to_geographic(x_min, x_max, y_min, y_max, out_gjson_path)



if __name__ == "__main__":

    #cc.file_query('/Users/roxana/OPEN/OPEN/boundary_info/s10.csv')
    cc.file_query('/home/extractor/hs_calib/s9_boundary.csv')
    input_paths = [
        # "scanDirectionIsPositive": "False" w-e
        #os.path.join(raw_root, "VNIR/2020-02-05/2020-02-05__11-39-58-853/a9dfc10e-b8c2-433e-9066-cd6830ce371f_raw"),
        #"scanDirectionIsPositive": "True" e-w
        # os.path.join(raw_root, "VNIR/2020-02-06/2020-02-06__12-39-41-588/17c5d1a7-ec83-4872-8fa9-df0b79cd8994_raw"),
        #"scanDirectionIsPositive": "False",
        os.path.join(raw_root, "VNIR/2019-06-09/2019-06-09__11-58-03-907/8fe14e6b-91a7-4637-aeff-e320e0bf5bff_raw"),
        #"scanDirectionIsPositive": "True"
        #os.path.join(raw_root, "VNIR/2019-06-09/2019-06-09__12-09-03-648/4133fb4e-75ed-49e7-8110-a8f0e7761cec_raw")
    ]

    for p in input_paths:
        process_VNIR(p, 1, out_root, env_log)
