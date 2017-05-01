#!/usr/bin/env python

"""
terra.hyperspectral.py

This extractor will trigger when a file is added to a dataset in Clowder.
It checks if all the required input files are present in the dataset while the
output file is not present. The output filename is always determined from the
filename of the `raw` file.
If the check is OK, it calls the `workerScript` defined in the config file to
create a netCDF output file and adds that to the same dataset.
"""

import os
import subprocess
import json
import logging
import tempfile

from pyclowder.extractors import Extractor
from pyclowder.utils import CheckMessage
import pyclowder.files
import pyclowder.datasets

class HyperspectralRaw2NetCDF(Extractor):
	def __init__(self):
		Extractor.__init__(self)

		# add any additional arguments to parser
		# self.parser.add_argument('--max', '-m', type=int, nargs='?', default=-1,
		#                          help='maximum number (default=-1)')
		self.parser.add_argument('--output', '-o', dest="output_dir", type=str, nargs='?',
								 default="/home/ubuntu/sites/ua-mac/Level_1/hyperspectral",
								 help="root directory where timestamp & output directories will be created")
		self.parser.add_argument('--script', dest="main_script", type=str, nargs='?',
								 default="hyperspectral_workflow.sh",
								 help="location of hyperspectral_workflow.sh file")

		# parse command line and load default logging configuration
		self.setup()

		# setup logging for the exctractor
		logging.getLogger('pyclowder').setLevel(logging.DEBUG)
		logging.getLogger('__main__').setLevel(logging.DEBUG)

		# assign other arguments
		self.output_dir = self.args.output_dir
		self.main_script = self.args.main_script

	def check_message(self, connector, host, secret_key, resource, parameters):
		if has_all_files(resource):
			if has_output_file(resource):
				logging.info('skipping dataset %s, output file already exists' % resource['id'])
				return CheckMessage.ignore
			else:
				# Check if we have necessary metadata, either as a .json file or attached to dataset
				found_md = False
				for f in resource['files']:
					if f['filename'] == 'metadata.json':
						found_md = True
				if not found_md:
					md = pyclowder.datasets.download_metadata(connector, host, secret_key,
															  resource['id'], self.extractor_info['name'])
					if len(md) > 0:
						for m in md:
							# Check if this extractor has already been processed
							if 'agent' in m and 'name' in m['agent']:
								if m['agent']['name'].find(self.extractor_info['name']) > -1:
									print("skipping dataset %s, already processed" % resource['id'])
									return CheckMessage.ignore
							if 'content' in m and 'lemnatec_measurement_metadata' in m['content']:
								found_md = True
				if found_md:
					return CheckMessage.download
		else:
			logging.info('skipping dataset %s, not all input files are ready' % resource['id'])
			return CheckMessage.ignore

	def process_message(self, connector, host, secret_key, resource, parameters):
		# Find input files in dataset
		target_files = {
			'raw': None,
			'raw.hdr': None,
			'image.jpg': None,
			'frameIndex.txt': None,
			'settings.txt': None,
			"_metadata.json": None
		}

		metafile = None
		ds_metafile = None
		last_path = None
		path_match = None
		tempdir = None
		symlinks = []
		for f in resource['local_paths']:
			for fileExt in target_files.keys():
				if f.endswith(fileExt) and fileExt != '_metadata.json':
					filedir = f.replace(os.path.basename(f), "")
					if not last_path:
						last_path = filedir
					else:
						if filedir != last_path:
							path_match = False
						last_path = filedir

					target_files[fileExt] = {
						'filename': os.path.basename(f),
						'path': f
					}
			if f.endswith('_metadata.json') and not f.endswith('/_metadata.json'):
				metafile = f
			if f.endswith('/_dataset_metadata.json'):
				ds_metafile = f

		# Identify md file either with other dataset files, or attached to Clowder dataset
		if metafile == None:
			if ds_metafile != None:
				# Found dataset metadata, so check for the .json file alongside other files
				logging.info("...checking for local metadata file alongside other files")
				ds_dir = target_files['raw']['path'].replace(os.path.basename(target_files['raw']['path']), "")
				ds_files = os.path.listdir(ds_dir)
				for ds_f in ds_files:
					if ds_f.endswith("_metadata.json"):
						target_files['_metadata.json']['path'] = os.path.join(ds_dir, ds_f)
			else:
				logging.error('could not locate metadata')
				return
		else:
			target_files['_metadata.json'] = {'filename': os.path.basename(metafile),
											  'path': metafile}

		# Create symlinks in one directory if inputs aren't in the same one
		if not path_match:
			tempdir = tempfile.mkdtemp()
			for f in target_files.keys():
				currf = target_files[f]
				if currf['filename'] == '_dataset_metadata.json':
					# Open the temporary file and change the JSON content format
					with open(currf['path'], 'r') as mdfile:
						jsondata = json.load(mdfile)
					for j in jsondata:
						if 'content' in j:
							if 'lemnatec_measurement_metadata' in j['content']:
								with open(currf['path'], 'w') as mdfile:
									json.dump(j['content'], mdfile)
					newf = os.path.join(tempdir, target_files['raw']['filename'].replace("_raw","")+'_metadata.json')
				else:
					newf = os.path.join(tempdir, currf['filename'])
				os.symlink(currf['path'], newf)
				symlinks.append(newf)

		# Prep output location
		try:
			date_portion = resource['dataset_info']['name'].split(' - ')[1].split('__')[0]
			timestamp_portion = resource['dataset_info']['name'].split(' - ')[1]
		except:
			date_portion = resource['dataset_info']['name']
			timestamp_portion = ""

		outFilePath = os.path.join(self.output_dir + ("_swir" if resource['dataset_info']['name'].find("SWIR") > -1
													  else ""),
								   date_portion,
								   timestamp_portion,
								   get_output_filename(target_files['raw']['filename']))
		out_dir = outFilePath.replace(os.path.basename(outFilePath), '')
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)

		# Invoke terraref.sh
		logging.debug('invoking terraref.sh to create: %s' % outFilePath)
		#returncode = subprocess.call(["bash", self.main_script, "-d", "2", "-I", inputDirectory, "-o", outFilePath])
		returncode = subprocess.call(["bash", self.main_script, "-d", "1", "-i", target_files['raw']['path'], "-o", outFilePath])

		# Verify outfile exists and upload to clowder
		logging.debug('done creating output file (%s)' % (returncode))
		if returncode != 0:
			logging.error('script encountered an error')
		if os.path.exists(outFilePath):
			if returncode == 0:
				logging.info('uploading %s' % outFilePath)
				pyclowder.files.upload_to_dataset(connector, host, secret_key, resource['id'], outFilePath)
		else:
			logging.error('no output file was produced')

		# Remove symlinks and temp directory
		for sym in symlinks:
			os.remove(sym)
		if tempdir:
			os.rmdir(tempdir)

# Find as many expected files as possible and return the set.
def get_all_files(resource):
	target_files = {
		'raw': None,
		'raw.hdr': None,
		'image.jpg': None,
		'frameIndex.txt': None,
		'settings.txt': None
	}

	if 'files' in resource:
		for fileItem in resource['files']:
			fileId   = fileItem['id']
			fileName = fileItem['filename']
			for fileExt in target_files.keys():
				if fileName.endswith(fileExt):
					target_files[fileExt] = {
						'id': fileId,
						'filename': fileName
					}

	return target_files

# Returns the output filename.
def get_output_filename(raw_filename):
	return '%s.nc' % raw_filename[:-len('raw')]

# Returns true if all expected files are found.
def has_all_files(resource):
	target_files = get_all_files(resource)

	allFilesFound = True
	for fileExt in target_files:
		if target_files[fileExt] == None:
			allFilesFound = False

	return allFilesFound

# Returns true if the output file is present.
def has_output_file(resource):
	if 'files' not in resource:
		return False
	if not has_all_files(resource):
		return False

	files = get_all_files(resource)
	outFilename = get_output_filename(files['raw']['filename'])
	outFileFound = False
	for fileItem in resource['files']:
		if outFilename == fileItem['filename']:
			outFileFound = True
			break

	return outFileFound

if __name__ == "__main__":
	extractor = HyperspectralRaw2NetCDF()
	extractor.start()
