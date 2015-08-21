#
#

 MAKE = make -j 4
#MAKE = make
#MAKE = make --debug

#MAKEOPT = -B

OPERATING_SYSTEM = $(shell uname)
ARCH = $(shell uname -m)

.PHONY: clean all debug optimize systeminfo guide
.PHONY: rebuild rebuild-debug rebuild-optimize ipo-opt

#all: debug optimize
all: debug ipo-opt

rebuild: rebuild-debug rebuild-optimize

std08: FCADD += -std08
std08: rebuild
std03: FCADD += -std03
std03: rebuild
std95: FCADD += -std95
std95: rebuild
std90: FCADD += -std90
std90: rebuild

warn-std08: FCADD += -std08 -warn stderrors
warn-std08: rebuild
warn-std03: FCADD += -std03 -warn stderrors
warn-std03: rebuild
warn-std95: FCADD += -std95 -warn stderrors
warn-std95: rebuild
warn-std90: FCADD += -std90 -warn stderrors
warn-std90: rebuild

debug:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=debug_recon" \
	                             "BUILDDIR=debug" \
				     "FCADD=$(FCADD)" \
				     debug
	@echo "-----------------------------------------------------"

optimize:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     optimize
	@echo "-----------------------------------------------------"

rebuild-debug:
	make -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=debug_recon" \
	                             "BUILDDIR=debug" \
				     "FCADD=$(FCADD)" \
				     debug
	@echo "-----------------------------------------------------"

rebuild-optimize:
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     optimize
	@echo "-----------------------------------------------------"

ipo-opt:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     ipo-opt

opt-report:
	$(MAKE) $(MAKEOPT) -B -f slv.mk "EXECUTABLE=recon" \
	                                "BUILDDIR=optimize" \
				        "FCADD=$(FCADD)" \
				        opt-report

mkl:
	$(MAKE) $(MAKEOPT) -B -f slv.mk "EXECUTABLE=recon" \
	                                "BUILDDIR=optimize" \
				        "FCADD=$(FCADD)" \
				        mkl

guide:
	$(MAKE) $(MAKEOPT) -B -f slv.mk "EXECUTABLE=recon" \
	                                "BUILDDIR=optimize" \
				        "FCADD=$(FCADD)" \
				        guide

mpe:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     "FC=mpefc -mpilog" \
				     optimize

mpe-check:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     "FC=mpefc -mpicheck -mpilog" \
				     optimize

systeminfo:
	@echo "Operating system = " $(OPERATING_SYSTEM)
	@echo "Architecture = " $(ARCH)

clean:
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean
	make -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" clean


