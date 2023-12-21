#!/bin/bash
#
# makestarmap.sh
#
# download and process a star map for use in Radiance with GenUtahSky
#

#
# Download the star map in EXR format
#

# get the 32k x 16k, 1.4GB version
#sudo dnf install curl
curl -LO https://svs.gsfc.nasa.gov/vis/a000000/a004800/a004851/starmap_2020_32k.exr

# 4k is way too low - the stars are unrealistically large, but it's only 34MB
#curl -LO https://svs.gsfc.nasa.gov/vis/a000000/a004800/a004851/starmap_2020_4k.exr

#
# Resize and convert to Radiance
#

# Imagemagick supports EXR
#sudo dnf install ImageMagick
# convert identify -list format | grep Radiance
#      HDR* HDR       rw+   Radiance RGBE image format
# convert identify -list format | grep EXR
#      EXR  EXR       rw-   High Dynamic-range (HDR) (OpenEXR 2.3.0)
convert starmap_2020_32k.exr -resize 50% starmap_2020_16k.hdr

# GDAL 3.1 supports EXR
#gdal_translate -ot Float -of GTiff -outsize 50% 50% starmap_2020_32k.exr starmap_2020_16k.tif
#ra_tiff -r starmap_2020_16k.tif > starmap_2020_16k.hdr

# Gimp supports EXR
#   load in Gimp, change to 32-bit floats,
#   reduce to 50% (to make smaller stars) using linear,
#   export as tiff, then
#ra_tiff -r reduced.tif > starfield.hdr

#
# Finally scale the exposure down to suit GenUtahSky
#

pfilt -1 -e -10 starmap_2020_16k.hdr > starmap_2020_16k_scaled.hdr


# This generates an hdr image that can be mapped to the sky dome at night
# and used as a replacement to TychoSkymapII.t5_08192x04096.hdr in stardome.rad
#
# see with
#pfilt -1 -e +9 -x /5 -y /5 TychoSkymapII.t5_08192x04096.hdr | ximage
#pfilt -1 -e +9 -x /10 -y /10 starmap_2020_16k_scaled.hdr | ximage
