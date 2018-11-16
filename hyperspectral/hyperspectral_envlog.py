import os
import subprocess

vnir_list = "/gpfs_scratch/vnir_netcdf/vnir_dates.txt"
env_log = "/projects/arpae/terraref/sites/ua-mac/Level_1/envlog_netcdf"
alt_out = "/gpfs_scratch/envlog_netcdf"

vnir_dates = []
with open(vnir_list, 'r') as inp:
    for date in inp:
        vnir_dates.append(date.rstrip())

# Prioritize dates w/ VNIR first, then the rest
for date in vnir_dates: #os.listdir(env_log)
    date_dir = os.path.join(env_log, date)

    if os.path.isdir(date_dir):
        allfiles = os.listdir(date_dir)
        ncfiles = [os.path.join(date_dir, f) for f in allfiles if f.endswith(".nc")]
        ncfiles.sort()

        out_file = "envlog_netcdf_L1_ua-mac_%s.nc" % date
        out_path = os.path.join(alt_out, date)
        if not os.path.isdir(out_path):
            os.makedirs(out_path)

        if not os.path.isfile(os.path.join(out_path, out_file)):
            print("concatenating %s files from %s" % (len(ncfiles), date))
            cmd = "ncrcat "+" ".join(ncfiles)+" "+os.path.join(out_path, out_file)
            subprocess.call([cmd], shell=True)
