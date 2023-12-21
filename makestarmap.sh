#!/bin/bash
#
# makestarmap.sh
#
# download and process a star map for use in Radiance with GenUtahSky
#

# get the 32k 1.4GB version
curl -LO https://svs.gsfc.nasa.gov/vis/a000000/a004800/a004851/starmap_2020_32k.exr

# 8k is way too low - the stars are unrealistically large, but it's only 34MB
#curl -LO https://svs.gsfc.nasa.gov/vis/a000000/a004800/a004851/starmap_2020_4k.exr

# Imagemagick supports EXR
# GDAL 3.1 supports EXR
# Gimp supports EXR
#   load in Gimp, change to 32-bit floats, reduce to 50% (to make smaller stars) using linear
#   export as tiff, then
#ra_tiff -r reduced.tif > starfield.hdr



