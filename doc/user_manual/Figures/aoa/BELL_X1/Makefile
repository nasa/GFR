FIGNAME := aoa_vectors

FIGURES := $(FIGNAME).png

.PHONY: all
all: $(FIGURES) clean

$(FIGURES): $(FIGNAME).tec \
            bell_x1_mesh.tec \
            bell_mx_mesh.lay \
            bell_mx_mesh.mcr
	tec360 -b --mesa bell_mx_mesh.mcr
	$(magickCMD) $(FIGNAME).png -shave 2x2 tmp.png
	/usr/bin/mv -f tmp.png $(FIGNAME).png
	/usr/bin/rm -f bell_x1_mesh.obj

$(FIGNAME).tec: $(FIGNAME).py
	python3 $(FIGNAME).py

bell_x1_mesh.tec: obj2tec.py bell_x1_mesh.obj
	python3 obj2tec.py

bell_x1_mesh.obj: bell_x1-low_resolution-obj.zip
	/usr/bin/unzip bell_x1-low_resolution-obj.zip
	/usr/bin/rm -f $(TMPFILES)

TMPFILES := bell_x1-full_resolution.mtl \
            bell_x1_diffuse.jpg \
	    bell_x1_normal_max.jpg \
	    bell_x1_occlusion.jpg 

magickCMD := /usr/bin/convert
#magickCMD := /usr/bin/magick

clean:
	/usr/bin/rm -f bell_x1_mesh.tec $(FIGNAME).tec batch.log

realclean: clean
	/usr/bin/rm -f $(FIGURES)
