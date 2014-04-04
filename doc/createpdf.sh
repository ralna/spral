#!/bin/bash

BASENAME=`echo "$1" | sed 's/\.tex//g'`
WRAPPERFILE=${BASENAME}_wrapper

echo "\\documentclass{spral}
\\begin{document}
\\include{$BASENAME}
\\end{document}
" > ${WRAPPERFILE}.tex

pdflatex ${WRAPPERFILE}.tex
pdflatex ${WRAPPERFILE}.tex
mv ${WRAPPERFILE}.pdf ${BASENAME}.pdf
rm ${WRAPPERFILE}.*
