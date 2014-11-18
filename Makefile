#
#

 MAKE = make -j 4
#MAKE = make
#MAKE = make --debug

#MAKEOPT = -B

OPERATING_SYSTEM = $(shell uname)
ARCH = $(shell uname -m)

all: both
#all: optimize
#all: debug

both:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" debug
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" optimize

debug:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" debug

optimize:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" optimize

guide:
	mkdir -p guide_files
	@echo $(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" guide

systeminfo:
	@echo "Operating system = " $(OPERATING_SYSTEM)
	@echo "Architecture = " $(ARCH)

clean:
	make -f slv.mk "EXECUTABLE=recon" clean
	make -f slv.mk "EXECUTABLE=debug_recon" clean


