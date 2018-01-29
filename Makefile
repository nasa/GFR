#
#

 MAKE = make -j 16
#MAKE = make
#MAKE = make --debug

#MAKEOPT = -B

OPERATING_SYSTEM = $(shell uname)
ARCH = $(shell uname -m)

#FCADD += -DDISABLE_PURE
#FCADD += -DPROFILE_ON
#FCADD += -DIGNORE_MPI_INFO

debug: FCADD += -DPROFILE_ON
rebuild-debug: FCADD += -DPROFILE_ON

.PHONY: all
all: debug_gfr
all: gfr
.PHONY: all-nodes
all-nodes: gfr_merope
all-nodes: gfr_westmere
all-nodes: gfr_sandybridge
all-nodes: gfr_ivybridge
all-nodes: gfr_broadwell
all-nodes: gfr_skylake

.PHONY: rebuild
rebuild: rebuild-debug
rebuild: rebuild-optimize

.PHONY: optimize
optimize: opt
.PHONY: rebuild-optimize
rebuild-optimize: rebuild-opt

debug_gfr: debug
gfr: optimize

gfr_skylake: skylake
gfr_broadwell: broadwell
gfr_haswell: haswell
gfr_ivybridge: ivybridge
gfr_sandybridge: sandybridge
gfr_westmere: westmere
gfr_merope: merope



# SKYLAKE

.PHONY: skylake
skylake: BDIR := optimize_skylake
skylake: EXEF := gfr_skylake
skylake: TRGT := skylake
skylake: build

.PHONY: rebuild-skylake
rebuild-skylake: BDIR := optimize_skylake
rebuild-skylake: EXEF := gfr_skylake
rebuild-skylake: TRGT := skylake
rebuild-skylake: clean-build
rebuild-skylake: build

.PHONY: sky
sky: skylake
.PHONY: rebuild-sky
rebuild-sky: rebuild-skylake

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
haswell: EXEF := gfr_haswell
haswell: TRGT := haswell
haswell: build

.PHONY: rebuild-haswell
rebuild-haswell: BDIR := optimize_haswell
rebuild-haswell: EXEF := gfr_haswell
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
ivybridge: EXEF := gfr_ivybridge
ivybridge: TRGT := ivybridge
ivybridge: build

.PHONY: rebuild-ivybridge
rebuild-ivybridge: BDIR := optimize_ivybridge
rebuild-ivybridge: EXEF := gfr_ivybridge
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
sandybridge: EXEF := gfr_sandybridge
sandybridge: TRGT := sandybridge
sandybridge: build

.PHONY: rebuild-sandybridge
rebuild-sandybridge: BDIR := optimize_sandybridge
rebuild-sandybridge: EXEF := gfr_sandybridge
rebuild-sandybridge: TRGT := sandybridge
rebuild-sandybridge: clean-build
rebuild-sandybridge: build

.PHONY: san
san: sandybridge
.PHONY: rebuild-san
rebuild-san: rebuild-sandybridge

# WESTMERE

.PHONY: westmere
westmere: BDIR := optimize_westmere
westmere: EXEF := gfr_westmere
westmere: TRGT := westmere
westmere: build

.PHONY: rebuild-westmere
rebuild-westmere: BDIR := optimize_westmere
rebuild-westmere: EXEF := gfr_westmere
rebuild-westmere: TRGT := westmere
rebuild-westmere: clean-build
rebuild-westmere: build

.PHONY: wes
wes: westmere
.PHONY: rebuild-wes
rebuild-wes: rebuild-westmere

# MEROPE

.PHONY: merope
merope: BDIR := optimize_merope
merope: EXEF := gfr_merope
merope: TRGT := merope
merope: FCADD += -DPROFILE_ON
merope: build

.PHONY: rebuild-merope
rebuild-merope: BDIR := optimize_merope
rebuild-merope: EXEF := gfr_merope
rebuild-merope: TRGT := merope
rebuild-merope: FCADD += -DPROFILE_ON
rebuild-merope: clean-build
rebuild-merope: build

.PHONY: mer
mer: merope
.PHONY: rebuild-mer
rebuild-mer: rebuild-merope


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

#gfr_haswell: $(call build-function, optimize_haswell, \
                       gfr_haswell, $(FCADD), haswell)
#rebuild-test: $(call rebuild-function, optimize_haswell, \
                      gfr_haswell, $(FCADD), haswell)
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
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=debug_gfr" \
	                             "BUILDDIR=debug" \
				     "FCADD=$(FCADD)" \
				     debug
	@echo "-----------------------------------------------------"

.PHONY: opt
opt:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=gfr" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     opt
	@echo "-----------------------------------------------------"

.PHONY: opt-host
opt-host:
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=gfr" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     opt-host
	@echo "-----------------------------------------------------"

.PHONY: rebuild-debug
rebuild-debug:
	make -f slv.mk "EXECUTABLE=debug_gfr" "BUILDDIR=debug" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=debug_gfr" \
	                             "BUILDDIR=debug" \
				     "FCADD=$(FCADD)" \
				     debug
	@echo "-----------------------------------------------------"

.PHONY: rebuild-opt
rebuild-opt:
	make -f slv.mk "EXECUTABLE=gfr" "BUILDDIR=optimize" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=gfr" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     opt
	@echo "-----------------------------------------------------"

.PHONY: rebuild-opt-host
rebuild-opt-host:
	make -f slv.mk "EXECUTABLE=gfr" "BUILDDIR=optimize" clean
	$(MAKE) $(MAKEOPT) -f slv.mk "EXECUTABLE=gfr" \
	                             "BUILDDIR=optimize" \
				     "FCADD=$(FCADD)" \
				     opt-host
	@echo "-----------------------------------------------------"

.PHONY: link-to-haswell
link-to-haswell:
	make -f slv.mk link-to-haswell

.PHONY: deplist
deplist:
	make -f slv.mk "EXECUTABLE=debug_gfr" "BUILDDIR=debug" clean
	$(MAKE) $(MAKEOPT) -B -f slv.mk "EXECUTABLE=debug_gfr" \
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
	make -f slv.mk "EXECUTABLE=debug_gfr" "BUILDDIR=debug" clean

.PHONY: clean-optimize
clean-optimize:
	make -f slv.mk "EXECUTABLE=gfr" "BUILDDIR=optimize" clean

.PHONY: clean-build
clean-build:
	make -f slv.mk "EXECUTABLE=$(EXEF)" "BUILDDIR=$(BDIR)" clean

.PHONY: clean
clean:
	make -f slv.mk "EXECUTABLE=debug_gfr" \
	               "BUILDDIR=debug" clean
	make -f slv.mk "EXECUTABLE=gfr" \
	               "BUILDDIR=optimize" clean
	make -f slv.mk "EXECUTABLE=gfr_merope" \
	               "BUILDDIR=optimize_merope" clean
	make -f slv.mk "EXECUTABLE=gfr_westmere" \
	               "BUILDDIR=optimize_westmere" clean
	make -f slv.mk "EXECUTABLE=gfr_sandybridge" \
	               "BUILDDIR=optimize_sandybridge" clean
	make -f slv.mk "EXECUTABLE=gfr_ivybridge" \
	               "BUILDDIR=optimize_ivybridge" clean
	make -f slv.mk "EXECUTABLE=gfr_haswell" \
	               "BUILDDIR=optimize_haswell" clean
	make -f slv.mk "EXECUTABLE=gfr_skylake" \
	               "BUILDDIR=optimize_skylake" clean
	if [ -L gfr_broadwell ] ; then /bin/rm -f gfr_broadwell ; fi

