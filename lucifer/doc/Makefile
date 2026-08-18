# Build the FEL physics manual: make          (needs a TeX installation with latexmk)
fel-physics.pdf: fel-physics.tex
	latexmk -pdf -quiet fel-physics.tex
clean:
	latexmk -C fel-physics.tex
