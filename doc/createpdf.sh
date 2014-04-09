#!/bin/bash

BASENAME=`echo "$1" | sed 's/\.tex//g'`
WRAPPERFILE=${BASENAME}_wrapper
WORKDIR=work.$$

mkdir $WORKDIR
pushd $WORKDIR

echo "\\documentclass{spral}
\\begin{document}
\\input{../$BASENAME}
\\end{document}
" > ${WRAPPERFILE}.tex

ln -s ../../spral.cls
ln -s ../../stfc.png
ln -s ../../epsrc.png
ln -s ../../../examples
pdflatex ${WRAPPERFILE}.tex
pdflatex ${WRAPPERFILE}.tex
mv ${WRAPPERFILE}.pdf ../${BASENAME}.pdf
popd
rm -r $WORKDIR
