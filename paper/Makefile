all : paper.pdf

paper.pdf : paper.md paper.bib example-modular-augur-workflows.pdf
	pandoc --citeproc --bibliography=paper.bib --template template.tex -s $< -o $@

example-modular-augur-workflows.pdf : example-modular-augur-workflows.dot
	dot -Tpdf $< > $@

clean :
	rm -f paper.pdf example-modular-augur-workflows.pdf
