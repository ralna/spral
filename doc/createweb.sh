#!/bin/bash

if test "x$DOCDIR" == "x"; then
   DOCDIR=/numerical/www/numerical-www/spral/doc/latest
fi

sphinx-build Fortran $DOCDIR/Fortran
sphinx-build C $DOCDIR/C
