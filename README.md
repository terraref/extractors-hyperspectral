# Hyperspectral extractors

This repository contains extractors that process data originating from:
- Hyperspec INSPECTOR SWIR camera
- Hyperspec INSPECTOR VNIR camera


### Hyperspectral extractor
This extractor processes HDF files into netCDF. 

_Input_

  - Evaluation is triggered whenever a file is added to a dataset
  - Checks whether the file is a _raw file
  
_Output_

  - The dataset containing the _raw file will get a corresponding .nc netCDF file
