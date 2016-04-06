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
.PHONY: variables-ivybridge

all: debug optimize
all-nodes: haswell ivybridge sandybridge westmere
rebuild: rebuild-debug rebuild-optimize
##clean-debug clean-optimize

optimize: O3-ipo-opt
rebuild-optimize: rebuild-O3-ipo-opt




#debug: BDIR := debug
#debug: EXEF := debug_recon
#debug: TRGT := debug
#debug: build

#rebuild-debug: BDIR := debug
#rebuild-debug: EXEF := debug_recon
#rebuild-debug: TRGT := debug
#rebuild-debug: build

haswell: BDIR := optimize_haswell
haswell: EXEF := recon_haswell
haswell: TRGT := haswell
haswell: build

ivybridge: BDIR := optimize_ivybridge
ivybridge: EXEF := recon_ivybridge
ivybridge: TRGT := ivybridge
ivybridge: build

sandybridge: BDIR := optimize_sandybridge
sandybridge: EXEF := recon_sandybridge
sandybridge: TRGT := sandybridge
sandybridge: build

westmere: BDIR := optimize_westmere
westmere: EXEF := recon_westmere
westmere: TRGT := westmere
westmere: build

rebuild-haswell: BDIR := optimize_haswell
rebuild-haswell: EXEF := recon_haswell
rebuild-haswell: TRGT := haswell
rebuild-haswell: clean
rebuild-haswell: build

rebuild-ivybridge: BDIR := optimize_ivybridge
rebuild-ivybridge: EXEF := recon_ivybridge
rebuild-ivybridge: TRGT := ivybridge
rebuild-ivybridge: clean
rebuild-ivybridge: build

rebuild-sandybridge: BDIR := optimize_sandybridge
rebuild-sandybridge: EXEF := recon_sandybridge
rebuild-sandybridge: TRGT := sandybridge
rebuild-sandybridge: clean
rebuild-sandybridge: build

rebuild-westmere: BDIR := optimize_westmere
rebuild-westmere: EXEF := recon_westmere
rebuild-westmere: TRGT := westmere
rebuild-westmere: clean
rebuild-westmere: build

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

build:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=$(EXEF)" \
	                             "BUILDDIR=$(BDIR)" \
				     "FCADD=$(FCADD)" \
				     $(TRGT)
	@echo "-----------------------------------------------------"

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

ivy_bridge: ivybridge
ivy-bridge: ivybridge

sandy_bridge: sandybridge
sandy-bridge: sandybridge

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

time-optimize:
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" time

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

clean:
	make -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" clean
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean
	make -f slv.mk "EXECUTABLE=recon_westmere" "BUILDDIR=optimize_westmere" clean
	make -f slv.mk "EXECUTABLE=recon_sandybridge" "BUILDDIR=optimize_sandybridge" clean
	make -f slv.mk "EXECUTABLE=recon_ivybridge" "BUILDDIR=optimize_ivybridge" clean
	make -f slv.mk "EXECUTABLE=recon_haswell" "BUILDDIR=optimize_haswell" clean

