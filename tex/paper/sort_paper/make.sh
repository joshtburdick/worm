#!/bin/bash
# Regenerates this.

tex2pdf figures.tex
unblurify < figures.pdf > figures1.pdf
mv figures1.pdf figures.pdf

tex2pdf supplementalFigures.tex
unblurify < supplementalFigures.pdf > supplementalFigures1.pdf
mv supplementalFigures1.pdf supplementalFigures.pdf

