# DictyPy

[![Build Status](https://travis-ci.org/kpj/DictyPy.svg?branch=master)](https://travis-ci.org/kpj/DictyPy)
[![Coverage Status](https://coveralls.io/repos/kpj/DictyPy/badge.svg?branch=master)](https://coveralls.io/r/kpj/DictyPy?branch=master)

Group all dicty genes into functional groups and compute respective codon usages.


## Usage

The main functionality can be found in `main.py`.
Furthermore, additional filters can be put into `filters.py`.


## Used Data

Dicty's primary CDS was acquired from [dictybase](http://dictybase.org/). The gene annotations were retrieved using [AmiGO 2](http://amigo.geneontology.org/amigo/landing).

## Used Modules

The external (non-standard) python modules used in this project are:
* [biopython](http://biopython.org/)
* [selenium](https://selenium-python.readthedocs.org/)
* [beautifulsoup](http://www.crummy.com/software/BeautifulSoup/bs4/doc/)
