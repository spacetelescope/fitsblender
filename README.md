# fitsblender

[![Build Status](https://dev.azure.com/spacetelescope/fitsblender/_apis/build/status/spacetelescope.fitsblender?branchName=master)](https://dev.azure.com/spacetelescope/fitsblender/_build/latest?definitionId=14&branchName=master)
[![codecov](https://codecov.io/gh/spacetelescope/fitsblender/branch/main/graph/badge.svg)](https://codecov.io/gh/spacetelescope/fitsblender)

FITSBLENDER
============

This package supports the creation of a combined header for a FITS file based on the contents of the headers
of a set of input FITS images.  A rules file defines what keywords will be present in the combined
output header as well as how the output value will be determined from the set of values from all 
the input image headers.  
