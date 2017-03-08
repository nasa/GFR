#
#

 MAKE = make -j 16
#MAKE = make
#MAKE = make --debug

#MAKEOPT = -B

OPERATING_SYSTEM = $(shell uname)
ARCH = $(shell uname -m)

#FCADD += -DPROFILE_ON


.PHONY: all
all: debug_recon
all: recon
.PHONY: all-nodes
all-nodes: recon_sandybridge
all-nodes: recon_ivybridge
all-nodes: recon_broadwell

.PHONY: rebuild
rebuild: rebuild-debug
rebuild: rebuild-optimize

.PHONY: optimize
optimize: O3-ipo-opt
.PHONY: rebuild-optimize
rebuild-optimize: rebuild-O3-ipo-opt

debug_recon: debug
recon: optimize

recon_broadwell: broadwell
recon_haswell: haswell
recon_ivybridge: ivybridge
recon_sandybridge: sandybridge

# BROADWELL (Essentially the same as Haswell)

.PHONY: broadwell
broadwell: haswell
broadwell: link-to-haswell

.PHONY: rebuild-broadwell
rebuild-broadwell: rebuild-haswell
rebuild-broadwell: link-to-haswell

.PHONY: bro
bro: broadwell
.PHONY: rebuild-bro
rebuild-bro: rebuild-broadwell

# HASWELL

.PHONY: haswell
haswell: BDIR := optimize_haswell
haswell: EXEF := recon_haswell
haswell: TRGT := haswell
haswell: build

.PHONY: rebuild-haswell
rebuild-haswell: BDIR := optimize_haswell
rebuild-haswell: EXEF := recon_haswell
rebuild-haswell: TRGT := haswell
rebuild-haswell: clean-build
rebuild-haswell: build

.PHONY: has
has: haswell
.PHONY: rebuild-has
rebuild-has: rebuild-haswell

# IVY BRIDGE

.PHONY: ivy_bridge
ivy_bridge: ivybridge
.PHONY: ivy-bridge
ivy-bridge: ivybridge

.PHONY: ivybridge
ivybridge: BDIR := optimize_ivybridge
ivybridge: EXEF := recon_ivybridge
ivybridge: TRGT := ivybridge
ivybridge: build

.PHONY: rebuild-ivybridge
rebuild-ivybridge: BDIR := optimize_ivybridge
rebuild-ivybridge: EXEF := recon_ivybridge
rebuild-ivybridge: TRGT := ivybridge
rebuild-ivybridge: clean-build
rebuild-ivybridge: build

.PHONY: ivy
ivy: ivybridge
.PHONY: rebuild-ivy
rebuild-ivy: rebuild-ivybridge

# SANDY BRIDGE

.PHONY: sandy_bridge
sandy_bridge: sandybridge
.PHONY: sandy-bridge
sandy-bridge: sandybridge

.PHONY: sandybridge
sandybridge: BDIR := optimize_sandybridge
sandybridge: EXEF := recon_sandybridge
sandybridge: TRGT := sandybridge
sandybridge: build

.PHONY: rebuild-sandybridge
rebuild-sandybridge: BDIR := optimize_sandybridge
rebuild-sandybridge: EXEF := recon_sandybridge
rebuild-sandybridge: TRGT := sandybridge
rebuild-sandybridge: clean-build
rebuild-sandybridge: build

.PHONY: san
san: sandybridge
.PHONY: rebuild-san
rebuild-san: rebuild-sandybridge


# CHECK AGAINST FORTRAN STANDARDS WITH ONLY COMPILER WARNINGS

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

# CHECK AGAINST FORTRAN STANDARDS TURNING WARNINGS TO COMPILE-TIME ERRORS

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

#recon_haswell: $(call build-function, optimize_haswell, recon_haswell, $(FCADD), haswell)
#rebuild-test: $(call rebuild-function, optimize_haswell, recon_haswell, $(FCADD), haswell)
define build-function
  $(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=$(1)" \
                               "BUILDDIR=$(2)" \
			       "FCADD=$(3)" \
			       $(4)
endef

define build-function
  $(call clean-function, $1, $2)
  $(call build-function, $1, $2, $3, $4)
endef

define clean-function
  make -f slv.mk "EXECUTABLE=$(1)" "BUILDDIR=$(2)" clean
endef

.PHONY: build
build:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=$(EXEF)" \
	                             "BUILDDIR=$(BDIR)" \
				     "FCADD=$(FCADD)" \
				     $(TRGT)
	@echo "-----------------------------------------------------"

.PHONY: debug
debug:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=debug_recon" \
	                             "BUILDDIR=debug" \
				     "FCADD=$(FCADD)" \
				     debug
	@echo "-----------------------------------------------------"

.PHONY: O2-ip-opt
O2-ip-opt:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     O2-ip-opt
	@echo "-----------------------------------------------------"

.PHONY: O3-ipo-opt
O3-ipo-opt:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     O3-ipo-opt
	@echo "-----------------------------------------------------"

.PHONY: rebuild-debug
rebuild-debug:
	make -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=debug_recon" \
	                             "BUILDDIR=debug" \
				     "FCADD=$(FCADD)" \
				     debug
	@echo "-----------------------------------------------------"

.PHONY: rebuild-O2-ip-opt
rebuild-O2-ip-opt:
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     O2-ip-opt
	@echo "-----------------------------------------------------"

.PHONY: rebuild-O3-ipo-opt
rebuild-O3-ipo-opt:
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     O3-ipo-opt
	@echo "-----------------------------------------------------"

.PHONY: time-optimize
time-optimize:
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" time

.PHOHY: link-to-haswell
link-to-haswell:
	/bin/ln -sf recon_haswell recon_broadwell
	/bin/ln -sf recon_haswell $(HOME)/bin/recon_broadwell

.PHONY: deplist
deplist:
	make -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" clean
	$(MAKE) $(MAKEOPT) -B -f slv.mk "EXECUTABLE=debug_recon" \
	                                "BUILDDIR=debug" \
				        "FCADD=$(FCADD)" \
				        deplist
	@echo "-----------------------------------------------------"

.PHONY: systeminfo
systeminfo:
	@echo "Operating system = " $(OPERATING_SYSTEM)
	@echo "Architecture = " $(ARCH)

.PHONY: clean-debug
clean-debug:
	make -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" clean

.PHONY: clean-optimize
clean-optimize:
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean

.PHONY: clean-build
clean-build:
	make -f slv.mk "EXECUTABLE=$(EXEF)" "BUILDDIR=$(BDIR)" clean

.PHONY: clean
clean:
	make -f slv.mk "EXECUTABLE=debug_recon" "BUILDDIR=debug" clean
	make -f slv.mk "EXECUTABLE=recon" "BUILDDIR=optimize" clean
	make -f slv.mk "EXECUTABLE=recon_sandybridge" "BUILDDIR=optimize_sandybridge" clean
	make -f slv.mk "EXECUTABLE=recon_ivybridge" "BUILDDIR=optimize_ivybridge" clean
	make -f slv.mk "EXECUTABLE=recon_haswell" "BUILDDIR=optimize_haswell" clean
