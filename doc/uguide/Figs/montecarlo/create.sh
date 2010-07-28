#!/bin/bash

fig2dev -L pstex $1.fig > $1.ps
fig2dev -L pstex_t -p $1.ps $1.fig > flowchart.tex
latex wrapper
dvips -t landscape -E wrapper.dvi -o $1.eps
ps2pdf @params.in $1.pdf