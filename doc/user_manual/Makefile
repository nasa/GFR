
#BuildMethod := rubber
 BuildMethod := pdfltx
#BuildMethod := latex
#BuildMethod := ltxmk

NAME := user_manual

MainFile := $(NAME).tex
CompiledFile := $(NAME).pdf

RUBBERPATH := ./rubber/bin/
RUBBEREXE := $(RUBBERPATH)rubber
RUBBERINFO := $(RUBBERPATH)rubber-info

BIBSRC = refs.bib

BIBINT = $(NAME).bbl
GLSINT = $(NAME).glo

BRFINT = $(NAME).brf

BIBINTS = $(BIBINT) $(NAME).blg
GLSINTS = $(GLSINT) $(NAME).acn $(NAME).idx

INTERMEDIATES = $(NAME).acr \
                $(NAME).alg \
                $(NAME).aux \
                $(NAME).glg \
                $(NAME).gls \
                $(NAME).ilg \
                $(NAME).ind \
                $(NAME).ist \
                $(NAME).log \
                $(NAME).lof \
                $(NAME).lot \
                $(NAME).lol \
                $(NAME).out \
                $(NAME).toc \
                $(NAME).rubbercache \
                $(NAME).listing \
                $(NAME).bcf \
                $(NAME).brf \
                $(NAME).run.xml \
                $(NAME).fdb_latexmk \
                $(NAME).fls \
                $(BIBINTS) \
                $(GLSINTS)

include latex/dep.mk

.PHONY: all
#all: print-file-lists
all: $(CompiledFile)
#all: rubber-build

.PHONY: rubber
rubber: $(FILES)
	$(RUBBEREXE) -f --pdf -s $(MainFile)

.PHONY: rubber-build
rubber-build: $(FILES)
	$(RUBBEREXE) -f --pdf -s $(MainFile)
	$(RUBBEREXE) -f --pdf -s $(MainFile)

$(CompiledFile): $(MainFile) $(FILES) $(BIBSRC)
	$(check-for-prebuild)
	$(build-pdf)

# Build standalone .pdf images from tikz .tex source files.

TIKZPDFS = tikz/EnergySpectra.pdf \
           tikz/Jump_Marker.pdf \
           tikz/couette_schematic.pdf \
           tikz/sem_partition_boundary.pdf \
           tikz/tikz_frames.pdf

$(TIKZPDFS): %.pdf: %.tex
	$(call build-tikz,$<)

define build-tikz
  @echo -e "\n Using pdflatex to pre-build PDF image from TikZ file $1...\n"
  pdflatex -output-directory tikz $1
endef

define check-for-prebuild
  $(if $(wildcard $(CompiledFile)),,$(call run-prebuild))
endef
define run-prebuild
  @echo -e "\n Running pre-build to fix rubber stopping prematurely...\n"
  $(if $(findstring rubber,$(BuildMethod)),$(call run-rubber),\
       $(if $(findstring pdfltx,$(BuildMethod)),$(call run-pdflatex),\
            $(if $(findstring latex,$(BuildMethod)),$(call run-latex),\
                 $(call run-pdflatex))))
endef

define run-rubber
  $(if $(wildcard $(RUBBEREXE)),\
    $(RUBBEREXE) -f --pdf -s $(MainFile),\
    $(call run-pdflatex))
endef
define run-pdflatex
  pdflatex $(NAME)
endef
define run-latex
  latex $(NAME)
endef

define build-pdf
  $(if $(findstring rubber,$(BuildMethod)),$(call rubber-check),\
       $(if $(findstring pdfltx,$(BuildMethod)),$(call use-pdflatex),\
            $(if $(findstring latex,$(BuildMethod)),$(call use-latex),\
                 $(if $(findstring ltxmk,$(BuildMethod)),$(call use-latexmk),\
                      $(call use-pdflatex)))))
endef

define rubber-check
  $(if $(wildcard $(RUBBEREXE)),\
       $(call use-rubber),\
       $(call use-pdflatex))
endef

define use-rubber
  @echo -e "\n Using rubber to build $(MainFile)...\n"
  $(comment pdflatex $(NAME))
  $(RUBBEREXE) -f --pdf -s $(MainFile)
endef

define use-pdflatex
  @echo -e "\n Using pdflatex to build $(MainFile)...\n"
  pdflatex       $(NAME)
  bibtex         $(NAME)
  $(comment makeglossaries $(NAME))
  pdflatex       $(NAME)
  pdflatex       $(NAME)
endef

define use-latex
  @echo -e "\n Using latex to build $(MainFile)...\n"
  latex $(NAME)
  bibtex $(NAME)
  latex $(NAME)
  latex $(NAME)
  dvips -t landscape $(NAME).dvi
  ps2pdf $(NAME).ps
endef

define use-latexmk
  @echo -e "\n Using latexmk to build $(MainFile)...\n"
  latexmk --pdf $(NAME)
endef

.PHONY: check
check:
	$(RUBBERINFO) --check $(NAME)

.PHONY: clean
clean :
	/bin/rm -f $(INTERMEDIATES) tikz/*.aux tikz/*.log

.PHONY: fullclean
fullclean : clean
	/bin/rm -f $(TIKZPDFS)

.PHONY: realclean
realclean : fullclean
	/bin/rm -f $(NAME).pdf

.PHONY: one
one: one-pdflatex

.PHONY: run
run: run-pdflatex

one-pdflatex: $(FILES)
	pdflatex  $(NAME)

run-pdflatex: one-pdflatex $(BIBINT)
	pdflatex  $(NAME)

$(BIBINT): $(BIBSRC)
	bibtex $(NAME)

$(GLSINT):
	makeglossaries $(NAME)

.PHONY: print-file-lists
print-file-lists:
	@echo -e ""
	@echo -e "NAME = $(NAME)\n"
	@echo -e "MainFile = $(MainFile)\n"
	@echo -e "FILES = $(FILES)\n"
	@echo -e "PhysModelsFiles = $(PhysModelsFiles)\n"
	@echo -e "BCsFiles = $(BCsFiles)\n"
	@echo -e "TimeSteppingFiles = $(TimeSteppingFiles)\n"
	@echo -e "NumericsFiles = $(NumericsFiles)\n"
	@echo -e "FilesFiles = $(FilesFiles)\n"
	@echo -e "NMLFiles = $(NMLFiles)\n"
	@echo -e "TestCasesFiles = $(TestCasesFiles)\n"

