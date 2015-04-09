# DictyPy
[![Build Status](https://travis-ci.org/kpj/DictyPy.svg?branch=master)](https://travis-ci.org/kpj/DictyPy)
[![Coverage Status](https://coveralls.io/repos/kpj/DictyPy/badge.svg?branch=master)](https://coveralls.io/r/kpj/DictyPy?branch=master)

Group all dicty genes into functional groups and compute respective codon usages.


## Usage
The program is started by calling `python main.py`. Blast-specific applications can be executed via `python blaster.py`.


## Notes for Blast
#### Setup
Setup a local blast instance by running `./blast_setup.sh`. Before initializing the local database, download the list of taxids-to-be-extracted via `python taxid_crawler.py`. Then simply execute `./blast_db_setup.sh` (this might take a while and requires loads of RAM).


## Notes for the Main Program
### Classifiers
Classifiers need to have a static variable `data_file` specifying the fasta file they are going to act on. They may furthermore have a static `preprocessing(genes)` function which acts on the given list of biopython's SeqRecord entries.
Their most important part is the `get_groupname(self, record)` method, which assigns a record to its logical group (as this function returns a list, it is possible to assign one entry to multiple groups). This group is then used in all later processing steps.

### Filters
All filters reside in `filters.py`. Each filter must have an `apply(self, record)` method, which takes one record (biopython's SeqRecord) as its argument.
This method returns `True` if the respective record should be kept and `False` otherwise.
A filter may have the static variable `post_annotation` (defaults to `False`), which declares that this filter should be applied after gene annotations took place.


## Used Data
Dicty's primary CDS was acquired from [dictybase](http://dictybase.org/). The gene annotations were retrieved using [AmiGO 2](http://amigo.geneontology.org/amigo/landing).

## Used Modules
The external (non-standard) python modules used in this project are:
* [biopython](http://biopython.org/)
* [selenium](https://selenium-python.readthedocs.org/)
* [beautifulsoup](http://www.crummy.com/software/BeautifulSoup/bs4/doc/)
