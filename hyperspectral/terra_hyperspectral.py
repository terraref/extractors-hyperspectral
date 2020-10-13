#!/usr/bin/env python

import os
import subprocess
from netCDF4 import Dataset

from pyclowder.utils import CheckMessage
from pyclowder.datasets import download_metadata, remove_metadata, upload_metadata
from pyclowder.files import upload_to_dataset, submit_extraction
from terrautils.metadata import get_extractor_metadata, get_terraref_metadata, get_season_and_experiment, \
	clean_metadata
from terrautils.extractors import TerrarefExtractor, is_latest_file, build_dataset_hierarchy_crawl, \
	contains_required_files, file_exists, load_json_file, check_file_in_dataset, build_metadata
from terrautils.betydb import submit_traits, add_arguments, get_site_boundaries

from hs_crop import process_VNIR


def add_local_arguments(parser):
	# add any additional arguments to parser
	# self.parser.add_argument('--max', '-m', type=int, nargs='?', default=-1,
	#                          help='maximum number (default=-1)')
	add_arguments(parser)

class HyperspectralRaw2NetCDF(TerrarefExtractor):
	def __init__(self):
		super(HyperspectralRaw2NetCDF, self).__init__()

		add_local_arguments(self.parser)

		# parse command line and load default logging configuration
		self.setup(sensor='vnir_netcdf')

		# assign other argumentse
		self.bety_url = self.args.bety_url
		self.bety_key = self.args.bety_key

	def check_message(self, connector, host, secret_key, resource, parameters):
		if "rulechecked" in parameters and parameters["rulechecked"]:
			return CheckMessage.download

		if not is_latest_file(resource):
			self.log_skip(resource, "not latest file")
			return CheckMessage.ignore

		if not contains_required_files(resource, ['raw', 'raw.hdr', 'image.jpg', 'frameIndex.txt', 'settings.txt']):
			self.log_skip(resource, "missing required files")
			return CheckMessage.ignore

		if resource['dataset_info']['name'].find("SWIR") > -1:
			sensor_fullname = 'swir_netcdf'
		else:
			sensor_fullname = 'vnir_netcdf'

		timestamp = resource['dataset_info']['name'].split(" - ")[1]
		md = download_metadata(connector, host, secret_key, resource['id'])
		if get_terraref_metadata(md):
			if get_extractor_metadata(md, self.extractor_info['name'], self.extractor_info['version']):
				# Make sure outputs properly exist
				out_nc = self.sensors.get_sensor_path(timestamp, sensor=sensor_fullname)
				if file_exists(out_nc):
					self.log_skip(resource, "metadata v%s and outputs already exist" % self.extractor_info['version'])
					return CheckMessage.ignore
			# Have TERRA-REF metadata, but not any from this extractor
			return CheckMessage.download
		else:
			self.log_skip(resource, "no terraref metadata found")
			# See if we can recover it from disk
			if sensor_fullname == 'vnir_netcdf':
				date = timestamp.split("__")[0]
				source_dir = "/home/extractor/sites/ua-mac/raw_data/VNIR/%s/%s/" % (date, timestamp)
				for f in os.listdir(source_dir):
					if f.endswith("_metadata.json"):
						self.log_info(resource, "updating metadata from %s" % f)
						raw_dsmd = load_json_file(os.path.join(source_dir, f))
						clean_md = clean_metadata(raw_dsmd, 'VNIR')
						complete_md = build_metadata(host, self.extractor_info, resource['id'], clean_md, 'dataset')
						remove_metadata(connector, host, secret_key, resource['id'])
						upload_metadata(connector, host, secret_key, resource['id'], complete_md)
						return CheckMessage.download
			return CheckMessage.ignore

	def process_message(self, connector, host, secret_key, resource, parameters):
		self.start_message(resource)

		# clean tmp directory from any potential failed previous runs
		flist = os.listdir("/tmp")
		for f in flist:
			try:
				os.remove(os.path.join("/tmp", f))
			except:
				pass

		timestamp = resource['dataset_info']['name'].split(" - ")[1]
		if resource['dataset_info']['name'].find("SWIR") > -1:
			sensor_rawname  = 'SWIR'
			sensor_fullname = 'swir_netcdf'
			soil_mask = None
		else:
			sensor_rawname  = 'VNIR'
			sensor_fullname = 'vnir_netcdf'
			# Check for corresponding soil mask to include in workflow.sh if available
			soil_mask = self.sensors.get_sensor_path(timestamp, sensor='vnir_soil_masks', opts=['soil_mask'])
		out_nc = self.sensors.create_sensor_path(timestamp, sensor=sensor_fullname)
		xps_file = self.sensors.get_sensor_path(timestamp, sensor=sensor_fullname, opts=['xps'])
		ind_file = self.sensors.get_sensor_path(timestamp, sensor=sensor_fullname, opts=['ind'])
		csv_file = self.sensors.get_sensor_path(timestamp, sensor=sensor_fullname.replace("_netcdf", "_traits"))

		raw_file, terra_md_full = None, None
		for fname in resource['local_paths']:
			if fname.endswith('_dataset_metadata.json'):
				all_dsmd = load_json_file(fname)
				terra_md_full = get_terraref_metadata(all_dsmd, sensor_rawname)
			elif fname.endswith('raw'):
				raw_file = fname
		if None in [raw_file, terra_md_full]:
			raise ValueError("could not locate all files & metadata in processing")

		# Fetch experiment name from terra metadata
		season_name, experiment_name, updated_experiment = get_season_and_experiment(timestamp, sensor_rawname, terra_md_full)
		if None in [season_name, experiment_name]:
			raise ValueError("season and experiment could not be determined")

		# Determine output directory
		print_name = self.sensors.get_display_name(sensor=sensor_fullname)
		self.log_info(resource, "Hierarchy: %s / %s / %s / %s / %s / %s / %s" % (season_name, experiment_name, print_name,
																				 timestamp[:4], timestamp[5:7], timestamp[8:10], timestamp))
		target_dsid = build_dataset_hierarchy_crawl(host, secret_key, self.clowder_user, self.clowder_pass, self.clowderspace,
													season_name, experiment_name, print_name,
													timestamp[:4], timestamp[5:7], timestamp[8:10],
													leaf_ds_name=self.sensors.get_display_name() + ' - ' + timestamp)
		uploaded_file_ids = []

		# Perform actual processing
		if (not file_exists(out_nc)) or self.overwrite:
			"""TODO: OLD AND NOT USED
			self.log_info(resource, 'invoking hyperspectral_workflow.sh to create: %s' % out_nc)
			if soil_mask and file_exists(soil_mask):
				# If soil mask exists, we can generate an _ind indices file
				returncode = subprocess.call(["bash", "hyperspectral_workflow.sh", "-d", "1", "-h",
										  "-m", soil_mask, "--output_xps_img", xps_file, "-i", raw_file, "-o", out_nc]) # disable --new_clb_mth
			else:
				# Otherwise we cannot, and need to trigger soilmask extractor and circle back later
				returncode = subprocess.call(["bash", "hyperspectral_workflow.sh", "-d", "1", "-h",
											  "--output_xps_img", xps_file, "-i", raw_file, "-o", out_nc]) # disable --new_clb_mth
			if returncode != 0:
				raise ValueError('script encountered an error')
			"""

			self.log_info(resource, 'invoking python calibration to create: %s' % out_nc)
			#create_empty_netCDF(raw_file, out_nc)
			#self.log_info(resource, 'applying calibration to: %s' % out_nc)
			#apply_calibration(raw_file, out_nc)

			# TODO: THIS IS NEW
			out_root = os.path.dirname(out_nc)
			env_root = "/home/extractor/sites/ua-mac/raw_data"
			process_VNIR(raw_file, 1, "/home/extractor/sites/ua-mac/Level_1_Plots/vnir_netcdf", env_root)

			self.log_info(resource, '...done' % raw_file)

			found_in_dest = check_file_in_dataset(connector, host, secret_key, target_dsid, out_nc, remove=self.overwrite)
			if not found_in_dest or self.overwrite:
				fileid = upload_to_dataset(connector, host, secret_key, target_dsid, out_nc)
				uploaded_file_ids.append(host + ("" if host.endswith("/") else "/") + "files/" + fileid)
			self.created += 1
			self.bytes += os.path.getsize(out_nc)

			# TODO: Soil mask still compatible in 3.0?
			#if not soil_mask:
			#	self.log_info(resource, "triggering soil mask extractor on %s" % fileid)
			#	submit_extraction(connector, host, secret_key, fileid, "terra.sunshade.soil_removal")

			# TODO: Send 3.0 output to BETYdb - needs more input
			"""
			# Send indices to betyDB
			if file_exists(ind_file):
				# TODO: Use ncks to trim ind_file to plots before this step
				plot_no = 'Full Field'
	
				with Dataset(ind_file, "r") as netCDF_handle:
					ndvi = netCDF_handle.get_variables_by_attributes(standard_name='normalized_difference_chlorophyll_index_750_705')
					NDVI705 = ndvi[0].getValue().ravel()[0]
	
					# TODO: Map the remaining ~50 variables in BETY to create indices file
					# TODO: In netCDF header,
	
					csv_header = 'local_datetime,NDVI705,access_level,species,site,' \
								 'citation_author,citation_year,citation_title,method'
					csv_vals = '%s,%s,2,Sorghum bicolor,%s,"Butowsky, Henry",2016,' \
							   'Maricopa Field Station Data and Metadata,Hyperspectral NDVI705 Indices' % (
									timestamp, NDVI705, plot_no)
					with open(csv_file, 'w') as c:
						c.write(csv_header+'\n'+csv_vals)
	
				# TODO: Send this CSV to betydb & geostreams extractors instead
				submit_traits(csv_file, bety_key=self.bety_key)
			"""

		self.end_message(resource)

if __name__ == "__main__":
	extractor = HyperspectralRaw2NetCDF()
	extractor.start()
