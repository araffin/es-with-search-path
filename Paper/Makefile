NAME=paper

all: *.tex
	@echo "make clean"
	@make clean
	@echo "pdflatex $(NAME).tex"
	@pdflatex $(NAME).tex
	# @echo "bibtex $(NAME).aux"
	# @bibtex $(NAME).aux
	@echo "pdflatex $(NAME).tex"
	@pdflatex $(NAME).tex
	@echo "pdflatex $(NAME).tex"
	@pdflatex $(NAME).tex
	@echo "make clean"
	@make clean

clean:
	rm -rf *.aux *.bbl *.blg *.glg *.glo *.gls *.idx  *.ilg *.ind *.ist *.log *.out *.toc *.backup *.xml *.bcf *.lot *.lof