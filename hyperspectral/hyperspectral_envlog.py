import os
import subprocess

env_log = "/projects/arpae/terraref/sites/ua-mac/Level_1/envlog_netcdf"
alt_out = "/gpfs_scratch/envlog_netcdf"

for date in os.listdir(env_log):
    date_dir = os.path.join(env_log, date)

    allfiles = os.listdir(date_dir)
    ncfiles = [f for f in allfiles if f.endswith(".nc")]
    ncfiles.sort()

    out_file = "envlog_netcdf_L1_ua-mac_%s.nc" % date
    out_path = os.path.join(alt_out, date)
    if not os.path.isdir(out_path):
        os.makedirs(out_path)

    print("concatenating %s files from %s" % (len(ncfiles), date))
    cmd = "ncrcat "+" ".join(ncfiles)+" "+os.path.join(out_path, out_file)
    subprocess.call([cmd], shell=True)
