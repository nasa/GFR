#
#

 MAKE = make -j 16
#MAKE = make -j 4
#MAKE = make
#MAKE = make --debug

#MAKEOPT = -B

OPERATING_SYSTEM = $(shell uname)
ARCH = $(shell uname -m)

FCADD += -DDISABLE_PURE
#FCADD += -DUSE_FACE_ARRAYS
FCADD += -DPROFILE_ON
#FCADD += -DIGNORE_MPI_INFO

#debug: FCADD += -DPROFILE_ON
#rebuild-debug: FCADD += -DPROFILE_ON

.PHONY: all
all: debug_recon
all: recon
.PHONY: all-nodes
all-nodes: recon_merope
all-nodes: recon_westmere
all-nodes: recon_sandybridge
all-nodes: recon_ivybridge
all-nodes: recon_broadwell
all-nodes: recon_skylake

.PHONY: rebuild
rebuild: rebuild-debug
rebuild: rebuild-optimize

.PHONY: optimize
optimize: O3-ipo-opt
.PHONY: rebuild-optimize
rebuild-optimize: rebuild-O3-ipo-opt

debug_recon: debug
recon: optimize

recon_skylake: skylake
recon_broadwell: broadwell
recon_haswell: haswell
recon_ivybridge: ivybridge
recon_sandybridge: sandybridge
recon_westmere: westmere
recon_merope: merope
recon_gpu_mm: gpu-mm
recon_gpu: gpu

# SKYLAKE

.PHONY: skylake
skylake: BDIR := optimize_skylake
skylake: EXEF := recon_skylake
skylake: TRGT := skylake
skylake: build

.PHONY: rebuild-skylake
rebuild-skylake: BDIR := optimize_skylake
rebuild-skylake: EXEF := recon_skylake
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

# WESTMERE

.PHONY: westmere
westmere: BDIR := optimize_westmere
westmere: EXEF := recon_westmere
westmere: TRGT := westmere
westmere: build

.PHONY: rebuild-westmere
rebuild-westmere: BDIR := optimize_westmere
rebuild-westmere: EXEF := recon_westmere
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
merope: EXEF := recon_merope
merope: TRGT := merope
#merope: FCADD += -pg
merope: build

.PHONY: rebuild-merope
rebuild-merope: BDIR := optimize_merope
rebuild-merope: EXEF := recon_merope
rebuild-merope: TRGT := merope
rebuild-merope: clean-build
rebuild-merope: build

.PHONY: mer
mer: merope
.PHONY: rebuild-mer
rebuild-mer: rebuild-merope

# GPU COMMON FLAGS

#FCADD += -acc -ta=tesla:managed,nollvm,ptxinfo,cc35 -Minfo=accel
#GPUFLAG = -acc -ta=tesla:nollvm,cc35,ptxinfo -Minfo=accel
GPUFLAG = -acc -ta=tesla:cc35,ptxinfo -Minfo=accel
#GPUADD += -Minfo=all
GPUADD += -Minfo=ccff,ftn

# GPU WITH MANAGED MEMORY

.PHONY: gpu-mm
gpu-mm: BDIR := gpu_mm
gpu-mm: EXEF := recon_gpu_mm
gpu-mm: TRGT := gpu-mm
gpu-mm: GPUADD += -ta=tesla:managed
#gpu-mm: GPUADD += -ta=tesla:fma
#gpu-mm: GPUADD += -ta=tesla:keepgpu
#gpu-mm: GPUADD += -ta=tesla:deepcopy
#gpu-mm: GPUADD += -ta=tesla:maxregcount:128
#gpu-mm: GPUADD += -Mnoidiom
#gpu-mm: GPUADD += -pg
#gpu-mm: GPUADD += -g
#gpu-mm: GPUADD += -Minfo=intensity
#gpu-mm: GPUADD += -Mcuda=lineinfo
#gpu-mm: GPUADD += -Mprof=ccff,dwarf
#gpu-mm: GPUADD += -Minstrument=functions
#gpu-mm: GPUADD += -Kieee
gpu-mm: FCADD += $(GPUFLAG) $(GPUADD)
gpu-mm: FCADD += -DDISABLE_FOR_GPU
gpu-mm: build

.PHONY: rebuild-gpu-mm
rebuild-gpu-mm: BDIR := gpu_mm
rebuild-gpu-mm: EXEF := recon_gpu_mm
rebuild-gpu-mm: TRGT := gpu-mm
rebuild-gpu-mm: clean-build
rebuild-gpu-mm: build

# GPU

.PHONY: gpu
gpu: BDIR := gpu
gpu: EXEF := recon_gpu
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
rebuild-gpu: EXEF := recon_gpu
rebuild-gpu: TRGT := gpu
rebuild-gpu: clean-build
rebuild-gpu: build


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

#recon_haswell: $(call build-function, optimize_haswell, \
                       recon_haswell, $(FCADD), haswell)
#rebuild-test: $(call rebuild-function, optimize_haswell, \
                      recon_haswell, $(FCADD), haswell)
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

.PHONY: link-to-haswell
link-to-haswell:
	make -f slv.mk link-to-haswell

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
	make -f slv.mk "EXECUTABLE=debug_recon" \
	               "BUILDDIR=debug" clean
	make -f slv.mk "EXECUTABLE=recon" \
	               "BUILDDIR=optimize" clean
	make -f slv.mk "EXECUTABLE=recon_merope" \
	               "BUILDDIR=optimize_merope" clean
	make -f slv.mk "EXECUTABLE=recon_westmere" \
	               "BUILDDIR=optimize_westmere" clean
	make -f slv.mk "EXECUTABLE=recon_sandybridge" \
	               "BUILDDIR=optimize_sandybridge" clean
	make -f slv.mk "EXECUTABLE=recon_ivybridge" \
	               "BUILDDIR=optimize_ivybridge" clean
	make -f slv.mk "EXECUTABLE=recon_haswell" \
	               "BUILDDIR=optimize_haswell" clean
	make -f slv.mk "EXECUTABLE=recon_skylake" \
	               "BUILDDIR=optimize_skylake" clean
	make -f slv.mk "EXECUTABLE=recon_gpu_mm" \
	               "BUILDDIR=gpu_mm" clean
	make -f slv.mk "EXECUTABLE=recon_gpu" \
	               "BUILDDIR=gpu" clean
	if [ -L recon_broadwell ] ; then /bin/rm -f recon_broadwell ; fi

