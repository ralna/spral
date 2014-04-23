#!/bin/bash

DOCDIR=/numerical/www/numerical-www/spral/doc

# Fortran
DESTDIR=$DOCDIR/Fortran/
WORKDIR=workF.$$
mkdir $WORKDIR
pushd $WORKDIR
ln -s ../Fortran/*.tex .
ln -s ../spral.tex ../install.tex ../spralweb.cls ../spralweb.cfg .
ln -s ../../examples
ln -s ../../LICENCE
htlatex spral.tex "spralweb.cfg, charset=utf-8" " -cunihtf -utf8" "-d$DESTDIR"
popd # WORKDIR
pushd $DESTDIR
ln -s spral.html index.html
popd # DESTDIR
rm -r $WORKDIR

# C
DESTDIR=$DOCDIR/C/
WORKDIR=workC.$$
mkdir $WORKDIR
pushd $WORKDIR
ln -s ../C/*.tex .
ln -s ../spral.tex ../install.tex ../spralweb.cls ../spralweb.cfg .
ln -s ../../examples
ln -s ../../LICENCE
htlatex spral.tex "spralweb.cfg, charset=utf-8" " -cunihtf -utf8" "-d$DESTDIR"
popd # WORKDIR
cp spral.custom.css spral.js $DESTDIR
cp stfc.png epsrc.png arrow_down16.png arrow_right16.png $DESTDIR
pushd $DESTDIR
ln -s spral.html index.html
popd # DESTDIR
rm -r $WORKDIR

# Shared
cp spral.custom.css spral.js $DOCDIR
cp stfc.png epsrc.png arrow_down16.png arrow_right16.png $DOCDIR
