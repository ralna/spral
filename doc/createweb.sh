#!/bin/bash

htlatex spral.tex "spralweb.cfg, charset=utf-8" " -cunihtf -utf8" "-d/numerical/www/numerical-www/spral/doc/"
rm *.html *.4ct *.4tc *.aux *.css *.dvi *.idv *.lg *.log *.tmp *.xref
