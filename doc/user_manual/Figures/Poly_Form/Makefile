FIGNAME := Poly_Form

FIGURES := $(FIGNAME)_Nodal.png $(FIGNAME)_Modal.png

.PHONY: all
all: $(FIGURES) clean

$(FIGURES): $(FIGNAME).tec \
            $(FIGNAME).exe \
            $(FIGNAME).lay \
            $(FIGNAME).mcr
	tec360 -b --mesa $(FIGNAME).mcr
	$(magickCMD) $(FIGNAME)_Nodal.png -shave 2x2 tmp.png
	/usr/bin/mv -f tmp.png $(FIGNAME)_Nodal.png
	$(magickCMD) $(FIGNAME)_Modal.png -shave 2x2 tmp.png
	/usr/bin/mv -f tmp.png $(FIGNAME)_Modal.png

$(FIGNAME).exe: $(FIGNAME).f90
	$(FC) $(FCFLAG) $^ -o $@

$(FIGNAME).tec: $(FIGNAME).exe
	echo "Running $(FIGNAME).exe"
	$(FIGNAME).exe

TMPFILES := $(FIGNAME).exe \
            batch.log

FC := gfortran
FCFLAG := 

magickCMD := /usr/bin/convert
#magickCMD := /usr/bin/magick

clean:
	/usr/bin/rm -f $(TMPFILES)

realclean: clean
	/usr/bin/rm -f $(FIGURES) $(FIGNAME).tec
