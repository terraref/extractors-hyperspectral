import os
import subprocess
import argparse
from terrautils.sensors import Sensors

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--input', help="input CSV file")
args = parser.parse_args()


sensor = Sensors('/projects/arpae/terraref/sites', 'ua-mac', 'vnir_netcdf')
vnir_raw = "/projects/arpae/terraref/sites/ua-mac/raw_data/VNIR"
vnir_lv1 = "/projects/arpae/terraref/sites/ua-mac/Level_1/vnir_netcdf"
alt_out = "/gpfs_scratch/vnir_netcdf"

with open(args.input, 'r') as inp:
    #for date in os.listdir(vnir_raw):
    for date in inp:
        date = date.rstrip()
        date_dir = os.path.join(vnir_raw, date)
        for timestamp in os.listdir(date_dir):
            # e.g. .../ua-mac/raw_data/VNIR/2016-04-05/2016-04-05__00-27-08-740/
            ds_dir = os.path.join(date_dir, timestamp)

            for hyperfile in os.listdir(ds_dir):
                if hyperfile.endswith("_raw"):
                    # Found a raw file - do we already have inputs?
                    lv1_out = sensor.get_sensor_path(timestamp)
                    if os.path.isdir(os.path.dirname(lv1_out)):
                        if not os.path.isfile(lv1_out):
                            # Invoke terraref.sh
                            script_path = "/projects/arpae/terraref/shared/extractors/extractors-hyperspectral/hyperspectral/hyperspectral_workflow.sh"
                            print(" ".join(["bash", script_path, "-d", "1", "-h", "--new_clb_mth",
                                   "-i", os.path.join(ds_dir, hyperfile),
                                   "-o", os.path.join(alt_out, date, timestamp, os.path.basename(lv1_out))]))
                            subprocess.call(["bash", script_path, "-d", "1", "-h", "--new_clb_mth",
                                                            "-i", os.path.join(ds_dir, hyperfile),
                                                            "-o", os.path.join(alt_out, date, timestamp, os.path.basename(lv1_out))])