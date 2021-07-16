#!/bin/sh

xelatex paper.tex
bibtex paper.aux
xelatex paper.tex
xelatex paper.tex
