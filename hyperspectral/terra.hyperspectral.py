#!/usr/bin/env python

"""
terra.hyperspectral.py

This extractor will trigger when a file is added to a dataset in Clowder.
It checks if all the required input files are present in the dataset while the
output file is not present. The output filename is always determined from the
filename of the `_raw` file.
If the check is OK, it calls the `workerScript` defined in the config file to
create a netCDF output file and adds that to the same dataset.
"""

import os
import subprocess
import json
import logging
from config import *
import pyclowder.extractors as extractors


def main():
	global extractorName, messageType, rabbitmqExchange, rabbitmqURL, registrationEndpoints, mountedPaths

	#set logging
	logging.basicConfig(format='%(levelname)-7s : %(name)s -  %(message)s', level=logging.WARN)
	logging.getLogger('pyclowder.extractors').setLevel(logging.INFO)
	logger = logging.getLogger('extractor')
	logger.setLevel(logging.DEBUG)

	print("main")
	print(rabbitmqURL)
	print(rabbitmqExchange)

	# setup
	extractors.setup(extractorName=extractorName,
					 messageType=messageType,
					 rabbitmqURL=rabbitmqURL,
					 rabbitmqExchange=rabbitmqExchange,
					 mountedPaths=mountedPaths)

	# register extractor info
	extractors.register_extractor(registrationEndpoints)

	#connect to rabbitmq
	extractors.connect_message_bus(extractorName=extractorName,
								   messageType=messageType,
								   processFileFunction=process_dataset,
								   checkMessageFunction=check_message,
								   rabbitmqExchange=rabbitmqExchange,
								   rabbitmqURL=rabbitmqURL)

def check_message(parameters):
	# Check for expected input files before beginning processing
	if has_all_files(parameters):
		if has_output_file(parameters):
			print 'skipping, output file already exists'
			return False
		else:
			# Check if we have necessary metadata, either as a .json file or attached to dataset
			found_md = False
			for f in parameters['filelist']:
				if f['filename'] == 'metadata.json':
					found_md = True
			if not found_md:
				md = extractors.download_dataset_metadata_jsonld(parameters['host'], parameters['secretKey'], parameters['datasetId'], extractorName)
				if len(md) > 0:
					for m in md:
						# Check if this extractor has already been processed
						if 'agent' in m and 'name' in m['agent']:
							if m['agent']['name'].find(extractorName) > -1:
								print("skipping dataset %s, already processed" % parameters['datasetId'])
								return False
						if 'content' in m and 'lemnatec_measurement_metadata' in m['content']:
							found_md = True
			if found_md:
				return True
	else:
		print 'skipping, not all input files are ready'
		return False

# ----------------------------------------------------------------------
# Process the dataset message and upload the results
def process_dataset(parameters):
	global extractorName, workerScript, inputDirectory, outputDirectory, requiredInputFiles

	# Find input files in dataset
	files = dict()
	metafile = None
	ds_metafile = None
	distinctPaths = []
	for f in parameters['files']:
		for fileExt in requiredInputFiles:
			if f.endswith(fileExt):
				files[fileExt] = {
					'filename': os.path.basename(f),
					'path': f
				}
				distinctPaths.append(f.replace(os.path.basename(f),''))
		if f.endswith('_metadata.json') and not f.endswith('/_metadata.json'):
			metafile = f
		if f.endswith('/_dataset_metadata.json'):
			ds_metafile = f

	# Identify md file either with other dataset files, or attached to Clowder dataset
	if metafile == None:
		if ds_metafile != None:
			files['_metadata.json'] = {'filename': os.path.basename(ds_metafile), 'path': ds_metafile}
			distinctPaths.append(metafile.replace(os.path.basename(metafile),''))
		else:
			print('could not locate metadata')
			return
	else:
		files['_metadata.json'] = {'filename': os.path.basename(metafile), 'path': metafile}
		distinctPaths.append(metafile.replace(os.path.basename(metafile),''))

	# Check if files are in same folder (i.e. is a local path) or not (i.e. is temp) - move together if not
	if len(distinctPaths) > 1:
		print("More than one directory found for files. Moving files to same directory for processing.")
		for fileExt in files:
			# Restore temp filenames to original - script requires specific name formatting so tmp names aren't suitable
			files[fileExt]['old_path'] = files[fileExt]['path']
			files[fileExt]['path'] = os.path.join(inputDirectory, files[fileExt]['filename'])
			os.rename(files[fileExt]['old_path'], files[fileExt]['path'])
			if fileExt == '_metadata.json' and files[fileExt]['filename'].endswith('_dataset_metadata.json'):
				# Convert _dataset_metadata autogenerated file into original metadata.json format
				rootName = files['_raw']['filename'].replace('_raw', '')
				newMetaName = files[fileExt]['path'].replace('_dataset_metadata.json', rootName+'_metadata.json')
				os.rename(files[fileExt]['path'], newMetaName)
				files[fileExt]['path'] = newMetaName
				with open(newMetaName, 'r') as md_file:
					md_contents = json.load(md_file)
				for md in md_contents:
					if 'content' in md and 'lemnatec_measurement_metadata' in md['content']:
						with open(newMetaName, 'w') as md_file:
							md_file.write(json.dumps(md['content']))
			print 'found %s file: %s' % (fileExt, files[fileExt]['path'])

	else:
		inputDirectory = distinctPaths[0]

	# Invoke terraref.sh
	outFilePath = os.path.join(outputDirectory,
							   parameters['datasetInfo']['name'].split(' - ')[1].split('__')[0],
							   parameters['datasetInfo']['name'].split(' - ')[1],
							   get_output_filename(files['_raw']['filename']))
	print 'invoking terraref.sh to create: %s' % outFilePath
	out_dir = outFilePath.replace(os.path.basename(outFilePath), '')
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	returncode = subprocess.call(["bash", workerScript, "-d", "2", "-I", inputDirectory, "-o", outFilePath])
	print 'done creating output file (%s)' % (returncode)

	if returncode != 0:
		print 'script encountered an error'

	# Verify outfile exists and upload to clowder
	if os.path.exists(outFilePath):
		if returncode == 0:
			print 'uploading output file...'
			extractors.upload_file_to_dataset(filepath=outFilePath, parameters=parameters)
		# Clean up the output file.
		os.remove(outFilePath)
	else:
		print 'no output file was produced'

	print 'cleaning up...'
	# Clean up the input files.
	for fileExt in files:
		os.remove(files[fileExt]['path'])

# ----------------------------------------------------------------------
# Find as many expected files as possible and return the set.
def get_all_files(parameters):
	global requiredInputFiles
	files = dict()
	for fileExt in requiredInputFiles:
		files[fileExt] = None

	if 'filelist' in parameters:
		for fileItem in parameters['filelist']:
			fileId   = fileItem['id']
			fileName = fileItem['filename']
			for fileExt in files:
				if fileName.endswith(fileExt):
					files[fileExt] = {
						'id': fileId,
						'filename': fileName
					}
	return files

# ----------------------------------------------------------------------
# Returns the output filename.
def get_output_filename(raw_filename):
	return '%s.nc' % raw_filename[:-len('_raw')]

# ----------------------------------------------------------------------
# Returns true if all expected files are found.
def has_all_files(parameters):
	files = get_all_files(parameters)

	allFilesFound = True
	for fileExt in files:
		if files[fileExt] == None:
			allFilesFound = False
	return allFilesFound

# ----------------------------------------------------------------------
# Returns true if the output file is present.
def has_output_file(parameters):
	if 'filelist' not in parameters:
		return False
	if not has_all_files(parameters):
		return False
	files = get_all_files(parameters)
	outFilename = get_output_filename(files['_raw']['filename'])
	outFileFound = False
	for fileItem in parameters['filelist']:
		if outFilename == fileItem['filename']:
			outFileFound = True
			break
	return outFileFound

if __name__ == "__main__":
	main()
