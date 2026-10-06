#!MC 1200

#
# Script to export the PNG image from a tecplot layout file.
#

$!PAGE NAME = 'Untitled'
$!PAGECONTROL CREATE

## Open tecplot layout file.
$!OPENLAYOUT  "bell_mx_mesh.lay"

## Set image export options.
$!PRINTSETUP PALETTE = COLOR
$!EXPORTSETUP EXPORTFORMAT = PNG
$!EXPORTSETUP IMAGEWIDTH   = 1504
$!EXPORTSETUP EXPORTFNAME = 'aoa_vectors.png'

## Export the image.
$!EXPORT
  EXPORTREGION = CURRENTFRAME

