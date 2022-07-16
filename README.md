# Extract_TROPOMI_Shapefiles

This repo is for code related to extracting data from TROPOMI NC files using shapefiles and/or clipping nc files using shapefiles.

## Notes

There is an error that may popup, but it doesn't halt the process and everything finishes find. It is likely related to garbage collector and can be safely ignored. Error in x$.self$finalize() : attempt to apply non-function

## Missing files

Just to double check that we got all of the files, as some files may have 0 soundings for our region of interest or we may have run out of memory, run check_missing.R to generate a list of potential missing files. Then run the main script again with that file list.
