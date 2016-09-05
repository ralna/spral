#!/bin/bash

DOCDIR=/numerical/www/numerical-www/spral/doc/sphinx

sphinx-build Fortran $DOCDIR/Fortran
sphinx-build C $DOCDIR/Fortran
