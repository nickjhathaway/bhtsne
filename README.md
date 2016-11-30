
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
./setup.py --libs armadillo:7.500.2 --outMakefile makefile-common.mk
```