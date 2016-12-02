
This software package contains a Barnes-Hut implementation of the t-SNE algorithm. The implementation is described in [this paper](http://lvdmaaten.github.io/publications/papers/JMLR_2014.pdf). This has now been edited by students at Umass Med in the BIB department


# Installation #


## Installing developmental version
```bash
git clone https://github.com/umass-bib/bhtsne.git
cd bhtsne 
git checkout develop
```

## Installing armadillo

Download and put makefile flags in makefile-common.mk

```bash
cd bhtsne 
./configure.py
./setup.py --compfile compfile.mk --outMakefile makefile-common.mk
```

Install with gcc compiler
```bash
cd bhtsne 
CC=gcc-5 CXX=g++-5 ./configure.py
./setup.py --compfile compfile.mk --outMakefile makefile-common.mk
```

