#!/bin/bash

DESTDIR=/numerical/www/numerical-www/spral/doc/

htlatex spral.tex "spralweb.cfg, charset=utf-8" " -cunihtf -utf8" "-d$DESTDIR"
rm *.html *.4ct *.4tc *.aux spral.css *.dvi *.idv *.lg *.log *.tmp *.xref
cp spral.custom.css spral.js $DESTDIR
cp stfc.png epsrc.png arrow_down16.png arrow_right16.png $DESTDIR
ln -s $DESTDIR/spral.html $DESTDIR/index.html
