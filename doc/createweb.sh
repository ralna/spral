#!/bin/bash

htlatex spral.tex "spralweb.cfg, 2, charset=utf-8, frames" " -cunihtf -utf8" "-d/numerical/www/numerical-www/spral/doc/"
rm *.html *.4ct *.4tc *.aux *.css *.dvi *.idv *.lg *.log *.tmp *.xref
