#
#

 MAKE = make -j 4
#MAKE = make
#MAKE = make --debug

#MAKEOPT = -B

OPERATING_SYSTEM = $(shell uname)
ARCH = $(shell uname -m)

.PHONY: all clean rebuild systeminfo deplist
.PHONY: debug      rebuild-debug     clean-debug
.PHONY: optimize   rebuild-optimize  clean-optimize
.PHONY: O2-ip-opt  rebuild-O2-ip-opt
.PHONY: O3-ipo-opt rebuild-O3-ipo-opt

all: debug optimize
rebuild: rebuild-debug rebuild-optimize
clean: clean-debug clean-optimize

optimize: O3-ipo-opt
rebuild-optimize: rebuild-O3-ipo-opt

.PHONY: std08
std08: FCADD += -std08
std08: rebuild
.PHONY: std03
std03: FCADD += -std03
std03: rebuild
.PHONY: std95
std95: FCADD += -std95
std95: rebuild
.PHONY: std90
std90: FCADD += -std90
std90: rebuild

.PHONY: warn-std08
warn-std08: FCADD += -warn stderrors
warn-std08: std08
.PHONY: warn-std03
warn-std03: FCADD += -warn stderrors
warn-std03: std03
.PHONY: warn-std95
warn-std95: FCADD += -warn stderrors
warn-std95: std95
.PHONY: warn-std90
warn-std90: FCADD += -warn stderrors
warn-std90: std90

#FCADD += -DTRACK_MEMORY_ALLOCATIONS
#FCADD += -DMEMORY_TEST_PAUSE

debug:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=debug_recon" \
	                             "BUILDDIR=debug" \
				     "FCADD=$(FCADD)" \
				     debug
	@echo "-----------------------------------------------------"

O2-ip-opt:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     O2-ip-opt
	@echo "-----------------------------------------------------"

O3-ipo-opt:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     O3-ipo-opt
	@echo "-----------------------------------------------------"

rebuild-debug:
	make -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=debug_recon" \
	                             "BUILDDIR=debug" \
				     "FCADD=$(FCADD)" \
				     debug
	@echo "-----------------------------------------------------"

rebuild-O2-ip-opt:
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     O2-ip-opt
	@echo "-----------------------------------------------------"

rebuild-O3-ipo-opt:
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     O3-ipo-opt
	@echo "-----------------------------------------------------"

deplist:
	make -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" clean
	$(MAKE) $(MAKEOPT) -B -f slv.mk "EXECUTABLE=debug_recon" \
	                                "BUILDDIR=debug" \
				        "FCADD=$(FCADD)" \
				        deplist
	@echo "-----------------------------------------------------"

systeminfo:
	@echo "Operating system = " $(OPERATING_SYSTEM)
	@echo "Architecture = " $(ARCH)

clean-debug:
	make -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" clean

clean-optimize:
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean


