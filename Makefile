#
#

 MAKE = make -j 16
#MAKE = make -j 4
#MAKE = make
#MAKE = make --debug

#MAKEOPT = -B

OPERATING_SYSTEM = $(shell uname)
ARCH = $(shell uname -m)

#FCADD += -DDISABLE_PURE
#FCADD += -DUSE_FACE_ARRAYS
#FCADD += -DPROFILE_ON
#FCADD += -DIGNORE_MPI_INFO

debug: FCADD += -DPROFILE_ON
rebuild-debug: FCADD += -DPROFILE_ON

FCADD += -DUSE_INTRINSIC_FINDLOC #-DPROFILE_ON

.PHONY: all
all: debug_gfr
all: gfr_westmere
#all: gfr
.PHONY: all-nodes
all-nodes: debug_gfr
all-nodes: gfr_westmere
all-nodes: gfr_skylake
all-nodes: gfr_rome
all-nodes: gfr

.PHONY: rebuild
rebuild: rebuild-debug
rebuild: rebuild-westmere

.PHONY: optimize
optimize: opt
.PHONY: rebuild-optimize
rebuild-optimize: rebuild-opt

debug_gfr: debug
gfr: optimize

gfr_rome: rome
gfr_skylake: skylake
gfr_haswell: haswell
gfr_westmere: westmere
gfr_gpu: gpu

# ROME

.PHONY: rome
rome: BDIR := optimize_rome
rome: EXEF := gfr_rome
rome: TRGT := rome
rome: build

.PHONY: rebuild-rome
rebuild-rome: BDIR := optimize_rome
rebuild-rome: EXEF := gfr_rome
rebuild-rome: TRGT := rome
rebuild-rome: clean-build
rebuild-rome: rome

.PHONY: clean-rome
clean-rome: BDIR := optimize_rome
clean-rome: EXEF := gfr_rome
clean-rome: clean-build

.PHONY: rom
rom: rome
.PHONY: rebuild-rom
rebuild-rom: rebuild-rome
.PHONY: clean-rom
clean-rom: clean-rome

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
rebuild-skylake: skylake

.PHONY: clean-skylake
clean-skylake: BDIR := optimize_skylake
clean-skylake: EXEF := gfr_skylake
clean-skylake: clean-build

.PHONY: sky
sky: skylake
.PHONY: rebuild-sky
rebuild-sky: rebuild-skylake
.PHONY: clean-sky
clean-sky: clean-skylake

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
rebuild-haswell: haswell

.PHONY: clean-haswell
clean-haswell: BDIR := optimize_haswell
clean-haswell: EXEF := gfr_haswell
clean-haswell: clean-build

.PHONY: has
has: haswell
.PHONY: rebuild-has
rebuild-has: rebuild-haswell
.PHONY: clean-has
clean-has: clean-haswell

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
rebuild-westmere: westmere

.PHONY: clean-westmere
clean-westmere: BDIR := optimize_westmere
clean-westmere: EXEF := gfr_westmere
clean-westmere: clean-build

.PHONY: wes
wes: westmere
.PHONY: rebuild-wes
rebuild-wes: rebuild-westmere
.PHONY: clean-wes
clean-wes: clean-westmere

# GPU COMMON FLAGS

#FCADD += -acc -ta=tesla:managed,nollvm,ptxinfo,cc35 -Minfo=accel
#GPUFLAG = -acc -ta=tesla:nollvm,cc35,ptxinfo -Minfo=accel
GPUFLAG = -acc -ta=tesla:cc35,ptxinfo -Minfo=accel
#GPUADD += -Minfo=all
GPUADD += -Minfo=ccff,ftn

# GPU

.PHONY: gpu
gpu: BDIR := gpu
gpu: EXEF := gfr_gpu
gpu: TRGT := gpu
#gpu: FCADD += -acc -ta=tesla:managed,nollvm,ptxinfo,cc35 -Minfo=accel
#gpu: GPUADD += -ta=tesla:nollvm
#gpu: GPUADD += -ta=tesla:fma
#gpu: GPUADD += -ta=tesla:keepgpu
#gpu: GPUADD += -ta=tesla:deepcopy
#gpu: GPUADD += -ta=tesla:maxregcount:128
#gpu: GPUADD += -Mnoidiom
#gpu: GPUADD += -pg
gpu: GPUADD += -g
#gpu: GPUADD += -Minfo=intensity
gpu: GPUADD += -ta=tesla:lineinfo
#gpu: GPUADD += -Mprof=ccff,dwarf
#gpu: GPUADD += -Minstrument=functions
#gpu: GPUADD += -Kieee
gpu: FCADD += $(GPUFLAG) $(GPUADD)
gpu: FCADD += -DDISABLE_FOR_GPU
gpu: build

.PHONY: rebuild-gpu
rebuild-gpu: BDIR := gpu
rebuild-gpu: EXEF := gfr_gpu
rebuild-gpu: TRGT := gpu
rebuild-gpu: clean-build
rebuild-gpu: gpu

.PHONY: clean-gpu
clean-gpu: BDIR := gpu
clean-gpu: EXEF := gfr_gpu
clean-gpu: clean-build


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

.PHONY: deplist
deplist:
	make -f slv.mk "EXECUTABLE=debug_gfr" "BUILDDIR=debug" clean
	make $(MAKEOPT) -B -f slv.mk "EXECUTABLE=debug_gfr" \
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
	make -f slv.mk "EXECUTABLE=gfr_westmere" \
	               "BUILDDIR=optimize_westmere" clean
	make -f slv.mk "EXECUTABLE=gfr_haswell" \
	               "BUILDDIR=optimize_haswell" clean
	make -f slv.mk "EXECUTABLE=gfr_skylake" \
	               "BUILDDIR=optimize_skylake" clean
	make -f slv.mk "EXECUTABLE=gfr_rome" \
	               "BUILDDIR=optimize_rome" clean
	make -f slv.mk "EXECUTABLE=gfr_gpu" \
	               "BUILDDIR=gpu" clean

