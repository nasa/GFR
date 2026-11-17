#
#

.NOTPARALLEL:

#MAKE := make

 MAKE_ARGS := -j 16
#MAKE_ARGS += --debug
#MAKE_ARGS += -B

#FCADD += -DDISABLE_PURE
#FCADD += -DUSE_FACE_ARRAYS
#FCADD += -DPROFILE_ON
#FCADD += -DIGNORE_MPI_INFO

#OPT_DEBUG_MODE := ON
#OPT_PROFILE_MODE := ON

.PHONY: all
all: debug opt

.PHONY: all-nodes
all-nodes: debug aitken athena opt

.PHONY: rebuild
rebuild: rebuild-debug rebuild-aitken

#######################
###  DEBUG TARGETS  ###
#######################
debug clean-debug rebuild-debug test-debug: BDIR := debug
debug clean-debug rebuild-debug test-debug: EXEF := debug_gfr
debug clean-debug rebuild-debug test-debug: TRGT := debug
debug clean-debug rebuild-debug test-debug: FCADD += -DDEBUG_ON -DPROFILE_ON

.PHONY: debug
debug:
	$(build-function)

.PHONY: clean-debug
clean-debug:
	$(clean-function)

.PHONY: rebuild-debug
#rebuild-debug: clean-debug debug
rebuild-debug:
	$(rebuild-function)

.PHONY: test-debug
test-debug:
	$(test-function)

#############################
###  HOST/NATIVE TARGETS  ###
#############################

.PHONY: opt-vars
opt clean-opt rebuild-opt test-opt: BDIR := optimize
opt clean-opt rebuild-opt test-opt: EXEF := gfr
opt clean-opt rebuild-opt test-opt: TRGT := opt
ifeq "$(OPT_DEBUG_MODE)" "ON"
opt clean-opt rebuild-opt test-opt: FCADD += -DDISABLE_PURE
opt clean-opt rebuild-opt test-opt: FCADD += -fpe0 -fpe-all=0 -fp-stack-check
endif
ifeq "$(OPT_PROFILE_MODE)" "ON"
opt clean-opt rebuild-opt test-opt: FCADD += -DPROFILE_ON
endif

.PHONY: opt
opt:
	$(build-function)

.PHONY: clean-opt
clean-opt:
	$(clean-function)

.PHONY: rebuild-opt
#rebuild-opt: clean-opt opt
rebuild-opt:
	$(rebuild-function)

.PHONY: test-opt
test-opt:
	$(test-function)

########################
###  AITKEN TARGETS  ###
########################

.PHONY: aitken-vars
aitken clean-aitken rebuild-aitken test-aitken: BDIR := optimize_aitken
aitken clean-aitken rebuild-aitken test-aitken: EXEF := gfr_aitken
aitken clean-aitken rebuild-aitken test-aitken: TRGT := aitken
ifeq "$(OPT_DEBUG_MODE)" "ON"
aitken clean-aitken rebuild-aitken test-aitken: FCADD += -DDISABLE_PURE
aitken clean-aitken rebuild-aitken test-aitken: FCADD += -fpe0 -fpe-all=0 -fp-stack-check
endif
ifeq "$(OPT_PROFILE_MODE)" "ON"
aitken clean-aitken rebuild-aitken test-aitken: FCADD += -DPROFILE_ON
endif

.PHONY: aitken
aitken:
	$(build-function)

.PHONY: clean-aitken
clean-aitken:
	$(clean-function)

.PHONY: rebuild-aitken
#rebuild-aitken: clean-aitken aitken
rebuild-aitken:
	$(rebuild-function)

.PHONY: test-aitken
test-aitken:
	$(test-function)

########################
###  ATHENA TARGETS  ###
########################

.PHONY: athena-vars
athena clean-athena rebuild-athena test-athena: BDIR := optimize_athena
athena clean-athena rebuild-athena test-athena: EXEF := gfr_athena
athena clean-athena rebuild-athena test-athena: TRGT := athena
ifeq "$(OPT_DEBUG_MODE)" "ON"
athena clean-athena rebuild-athena test-athena: FCADD += -DDISABLE_PURE
athena clean-athena rebuild-athena test-athena: FCADD += -fpe0 -fpe-all=0 -fp-stack-check
endif
ifeq "$(OPT_PROFILE_MODE)" "ON"
athena clean-athena rebuild-athena test-athena: FCADD += -DPROFILE_ON
endif

.PHONY: athena
athena:
	$(build-function)

.PHONY: clean-athena
clean-athena:
	$(clean-function)

.PHONY: rebuild-athena
#rebuild-athena: clean-athena athena
rebuild-athena:
	$(rebuild-function)

.PHONY: test-athena
test-athena:
	$(test-function)

#####################
###  GPU TARGETS  ###
#####################

.PHONY: gpu-vars
gpu clean-gpu rebuild-gpu test-gpu: BDIR := gpu
gpu clean-gpu rebuild-gpu test-gpu: EXEF := gfr_gpu
gpu clean-gpu rebuild-gpu test-gpu: TRGT := gpu
#gpu clean-gpu rebuild-gpu test-gpu: \
  GPUFLAG = -acc -ta=tesla:nollvm,cc35,ptxinfo -Minfo=accel
gpu clean-gpu rebuild-gpu test-gpu: GPUFLAG = -acc -ta=tesla:cc35,ptxinfo -Minfo=accel
#gpu clean-gpu rebuild-gpu test-gpu: \
  GPUADD += -acc -ta=tesla:managed,nollvm,ptxinfo,cc35 -Minfo=accel
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -Minfo=all
gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -Minfo=ccff,ftn -g -ta=tesla:lineinfo
#gpu clean-gpu rebuild-gpu test-gpu: \
  FCADD += -acc -ta=tesla:managed,nollvm,ptxinfo,cc35 -Minfo=accel
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -ta=tesla:nollvm
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -ta=tesla:fma
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -ta=tesla:keepgpu
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -ta=tesla:deepcopy
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -ta=tesla:maxregcount:128
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -Mnoidiom
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -pg
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -Minfo=intensity
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -Mprof=ccff,dwarf
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -Minstrument=functions
#gpu clean-gpu rebuild-gpu test-gpu: GPUADD += -Kieee
gpu clean-gpu rebuild-gpu test-gpu: FCADD += $(GPUFLAG) $(GPUADD) -DDISABLE_FOR_GPU
ifeq "$(OPT_DEBUG_MODE)" "ON"
gpu clean-gpu rebuild-gpu test-gpu: FCADD += -DDISABLE_PURE
gpu clean-gpu rebuild-gpu test-gpu: FCADD += -fpe0 -fpe-all=0 -fp-stack-check
endif
ifeq "$(OPT_PROFILE_MODE)" "ON"
gpu clean-gpu rebuild-gpu test-gpu: FCADD += -DPROFILE_ON
endif

.PHONY: gpu
gpu:
	$(build-function)

.PHONY: clean-gpu
clean-gpu:
	$(clean-function)

.PHONY: rebuild-gpu
#rebuild-gpu: clean-gpu gpu
rebuild-gpu:
	$(rebuild-function)

.PHONY: test-gpu
test-gpu:
	$(test-function)

###############################
###  Function definitions  ####
###############################

define clean-function
	make -f slv.mk "EXECUTABLE=$(EXEF)" "BUILDDIR=$(BDIR)" clean
endef

define build-function
	$(MAKE) $(MAKE_ARGS) -f slv.mk "EXECUTABLE=$(EXEF)" \
				       "BUILDDIR=$(BDIR)" \
				       "FCADD=$(FCADD)" \
				       $(TRGT)
	@echo "-----------------------------------------------------"
endef

define rebuild-function
	$(clean-function)
	$(build-function)
endef

define test-function
	$(MAKE) $(MAKE_ARGS) -f slv.mk "EXECUTABLE=$(EXEF)" \
				       "BUILDDIR=$(BDIR)" \
				       "FCADD=$(FCADD)" \
				       dump-variables
	@echo "-----------------------------------------------------"
endef

######################################
###  Old build and clean targets  ####
######################################

.PHONY: build
build:
	$(MAKE) $(MAKE_ARGS) -f slv.mk "EXECUTABLE=$(EXEF)" \
				       "BUILDDIR=$(BDIR)" \
				       "FCADD=$(FCADD)" \
				       $(TRGT)
	@echo "-----------------------------------------------------"

.PHONY: dump-variables
dump-variables:
	$(MAKE) $(MAKE_ARGS) -f slv.mk "EXECUTABLE=$(EXEF)" \
				       "BUILDDIR=$(BDIR)" \
				       "FCADD=$(FCADD)" \
				       dump-variables
	@echo "-----------------------------------------------------"

.PHONY: systeminfo
systeminfo:
	@echo "Operating system = " $(shell uname)
	@echo "Architecture = " $(shell uname -m)

.PHONY: clean-build
clean-build:
	make -f slv.mk "EXECUTABLE=$(EXEF)" "BUILDDIR=$(BDIR)" clean

.PHONY: clean
clean: clean-debug clean-opt clean-aitken clean-athena clean-gpu

#########################################
###  ALIASES FOR LONGER TARGET NAMES  ###
#########################################

.PHONY: debug_gfr gfr_aitken gfr_athena
debug_gfr: debug
gfr_aitken: aitken
gfr_athena: athena

#=================================================================

.PHONY: dbg rebuild-dbg clean-dbg test-dbg
dbg: debug
rebuild-dbg: rebuild-debug
clean-dbg: clean-debug
test-dbg: test-debug

#=================================================================

.PHONY: ait rebuild-ait clean-ait test-ait
.PHONY: rome rebuild-rome clean-rome test-rome
.PHONY: rom rebuild-rom clean-rom test-rom
ait rome rom: aitken
rebuild-ait rebuild-rome rebuild-rom: rebuild-aitken
clean-ait clean-rome clean-rom: clean-aitken
test-ait test-rome test-rom: test-aitken

#=================================================================

.PHONY: ath rebuild-ath clean-ath test-ath
.PHONY: turin rebuild-turin clean-turin test-turin
.PHONY: tur rebuild-tur clean-tur test-tur
ath turin tur: athena
rebuild-ath rebuild-turin rebuild-tur: rebuild-athena
clean-ath clean-turin clean-tur: clean-athena
test-ath test-turin test-tur: test-athena

#=================================================================

.PHONY: host rebuild-host clean-host test-host
.PHONY: optimize rebuild-optimize clean-optimize test-optimize
.PHONY: gfr rebuild-gfr clean-gfr test-gfr
host optimize gfr: opt
rebuild-host rebuild-optimize rebuild-gfr: rebuild-opt
clean-host clean-optimize clean-gfr: clean-opt
test-host test-optimize test-gfr: test-opt

#=================================================================
