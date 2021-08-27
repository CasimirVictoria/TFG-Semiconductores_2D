#!/bin/bash
jupytext ../Jupyter-Book/formulacio_matriu_dinamica.md -o formulacio_matriu_dinamica.py && 
mv  formulacio_matriu_dinamica.py formulacio_matriu_dinamica.sage &&
pdflatex --shell-escape TFG-Casimir.tex &&
bibtex TFG-Casimir &&
bibtex TFG-Casimir &&
sage TFG-Casimir.sagetex.sage &&
pdflatex --shell-escape TFG-Casimir.tex
