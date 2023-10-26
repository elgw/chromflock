# CHROMFLOCK

![example structure](doc/cf_000038.png)

## Introduction

**chromflock** is a software package to deconvolve bulk Hi-C data into
putative single-cell structures. The overall scheme is borrowed from
[PGS](https://www.github.com/alberlab/PGS) although it is not a direct
fork.

Chromflock will:
 - Work with haploid as well as diploid structures.
 - Can integrate GPSeq data for radial preferences.
 - Supports spherical as well as an ellipsoidal domain for the beads.

The documentation is not complete at the moment and it is suggested
that anyone interested in chromflock start by reading [Kalhor et al,
2012](https://doi.org/10.1038/nbt.2057), especially the supplementary
materials since much of the terminology used here can be traced back
to that paper.

## Installation
There are no pre built packages for chromflock so it has to be built
from source. If all dependencies are installed it should be simple as:

``` shell
make
./makedeb-ubuntu_2204.sh
sudo apt-get install ./chromflock_x.y.z_amd64.deb
```
more details can be found in [INSTALL.md](INSTALL.md).

## Usage
Chromflock requires at least two data inputs

1. A contact probability matrix, $`A`$, where $`A(i,j)`$ is the
probability that bead $`i`$ and $`j`$ is in contact. This is typically
constructed from a (bulk) Hi-C matrix.
2. A label vector, $`L`$ defining how many beads there are and what
   numerical label each bead has.

and produce a user-specified number of putative single-cell 3D
structures with a binary contact matrix for details, see
[USAGE.md](USAGE.md) as well as the man pages.

## References

chromflock was used in
 - [Girelli et al. 2020](https://www.nature.com/articles/s41587-020-0519-y).

The ideas behind it can be found in the following papers:
 - [Kalhor et al, 2012](https://doi.org/10.1038/nbt.2057)
 - [Tjong et al, 2016](http://dx.doi.org/10.1073/pnas.1512577113)
 - [Hua et al, 2018](http://dx.doi.org/10.1038/nprot.2018.008)
   [github](https://www.github.com/alberlab/PGS)

For random numbers chromflock uses
 -  [McFarland,
    2014](http://www.tandfonline.com/doi/abs/10.1080/00949655.2015.1060234)
    [code](https://github.com/cd-mcfarland/fast_prng).
