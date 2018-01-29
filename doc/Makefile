NAME=inputs


FILES = $(NAME).tex

INTERMEDIATES = $(NAME).aux \
                $(NAME).log \
                $(NAME).ps \
                $(NAME).dvi

.PHONY: all
all: use-rubber

.PHONY: use-latexmk
use-latexmk: $(FILES)
	latexmk --pdf $(NAME)

.PHONY: use-pdflatex
use-pdflatex: $(FILES)
	pdflatex $(NAME)

.PHONY: use-rubber
use-rubber: $(FILES)
	rubber -f --pdf -s $(NAME).tex

.PHONY: use-latex
use-latex: $(FILES)
	latex $(NAME)
	latex $(NAME)
	latex $(NAME)
	dvips -t legal -t landscape $(NAME).dvi
	ps2pdf $(NAME).ps

.PHONY: check
check:
	rubber-info --check $(NAME)


.PHONY: clean
clean :
	/bin/rm -f $(INTERMEDIATES)

.PHONY: realclean
realclean :
	/bin/rm -f $(INTERMEDIATES) $(NAME).pdf


