#
#

MAKE = make -j 4

all:
	$(MAKE) -B -f slv.mk "EXECUTABLE=debug_recon" debug
	$(MAKE) -B -f slv.mk "EXECUTABLE=recon" optimize

debug:
	$(MAKE) -B -f slv.mk "EXECUTABLE=debug_recon" debug

optimize:
	$(MAKE) -B -f slv.mk "EXECUTABLE=recon" optimize

clean:
	make -f slv.mk "EXECUTABLE=recon" clean
	make -f slv.mk "EXECUTABLE=debug_recon" clean


