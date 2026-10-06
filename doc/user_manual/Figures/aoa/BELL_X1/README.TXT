Purpose:

- This directory contains files for plotting the angle of attack vectors
  onto the Bell-X1 airplane using Tecplot.

- Downloaded bell_x1-low_resolution-obj.zip from the Smithsonian.
    https://3d.si.edu/object/3d/bell-x-1:6c69a6bb-55e6-4356-8725-120ff7f8d652


Reason for Keeping:

- These files are needed for creating an image for the GFR user manual.


Instructions:

```
Simply run the command 'make'


Below are the individual commands run by the Makefile
-----------------------------------------------------

#!/bin/csh

  unzip bell_x1-low_resolution-obj.zip

# Convert the geometry data to Tecplot.
  if ( ! -e bell_x1_mesh.tec ) then
    obj2tec.py
  endif

# Create Tecplot file with the aoa vectors.
  if ( ! -e aoa_vectors.tec ) then
    aoa_vectors.py
  endif

# Create Tecplot image.
  tec360
    - Load:  bell_mx_mesh.lay
    - Export PNG image:
      - Image width = 1504
      - Antialiasing = On
        - Supersample factor = 3
      - Filename = aoa_vectors.png

# Use Imagemagick to shave the image to remove border from antialiasing.
  magick aoa_vectors.png -shave 2x2  aoa_vectors_shaved.png
  mv -f aoa_vectors_shaved.png aoa_vectors.png

```

Misc:

- Contents of bell_x1-low_resolution-obj.zip

```
Archive:  bell_x1-low_resolution-obj.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
      211  02-21-2020 13:52   bell_x1-full_resolution.mtl
  7338420  02-21-2020 13:44   bell_x1_diffuse.jpg
 11172137  02-21-2020 13:52   bell_x1_mesh.obj
  3018618  02-21-2020 13:44   bell_x1_normal_max.jpg
  4179827  02-21-2020 13:44   bell_x1_occlusion.jpg
---------                     -------
 25709213                     5 files
```

- The intermediate files can be removed if necessary.

