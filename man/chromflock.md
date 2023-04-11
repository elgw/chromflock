% CHROMFLOCK(1) Version ABCDE | chromflock documentation
% Erik Wernersson
% 2023

# NAME
**chromflock** is a set of tools to deconvolve bulk Hi-C data into
individual contact indication matrices and generate 3D structures.

# SYNOPSIS
The following commands are available:

- **init** initialize a folder for chromflock, see separate section below.
- **hic2cpm** convert Hi-C data to a chromflock contact probability
  matrix, see separate section below.
- **sprite2cpm**
- **string2any** can be used to quickly create small examples of input
  to **aflock** and **mflock**
- **any2string** can be used to quickly inspect the input files to **chromflock**


# INIT

## EXAMPLE

``` shell
 $ mkdir 10000structures
 $ cd 10000structures
 $ # create chromflock_gen
 $ chromflock init
 $ # edit settings and point to correct A-file and possible G-file.
 $ nvim chromflock_gen
 # Generate a (linear) batch script with everything to be run
 $ ./chromflock_gen
 # Run everything:
 $ ./chromflock_run
```

## TIPS
Use screen or some other program to run chromflock in the background.
This allows you to log out from a (remote) computer while chromflock
is still running.

chromflock_run outputs some progress information to status.txt while
running.

# HIC2CPM
**hic2cpm** Can be used to convert a Hi-C matrix into a contact
probability matrix (CPM) suitable for chromflock.

In this context a HiC-matrix is Hi-C data mapped and binned to a
genome and saved as a dense matrix in 64-bit floating point format
(double). The values of the Hi-C matrix are either raw or normalized
counts and hence the range of the values does not matter to
**hic2cpm**.

Further downstream the CPM matrix can be read by **aflock**. The
elements of the CPM are interpreted so that if `CPM_{ij}=p`, the
probability that bin `i` and `j` are in contact is `p`. I.e. if
`p==0`, `i` and `j` should never be in contact in any structure. On
the other side if `p==1` `i` and `j` should be in contact in all
structures in the population.

## METHOD
 - An (optional) matrix balancing is applied to give equal number of
   contacts per bead/bin, i.e. a **Doubly Stochastic Matrix**. For
   this the Sinkhorn–Knopp Algorithm (_SKA_) is used.
 - The matrix is scaled to reach the specific maximum number of
  contacts per bead that was specified with **--nCont**.
 - The first off-diagonal is set to 1 for all bins within a
   chromosome. This is to maintain connection between adjacent beads
   on the same chromosome.
 - Chromosome Y (or whatever is labelled as 24) is removed unless
   **--y** is specified.

The number of contacts per bead (in average) is important due to
geometric constraints. For a dense sphere packing
[https://en.wikipedia.org/wiki/Sphere_packing] the number of neighbors
per bead can not exceed 12, and two are already in use for the
previous and next bead (unless the bead is at the beginning or end of
a chromosome). Since we don't expect perfect lattices in our
structures the number of contacts per bead should be strictly below 12.

Obviously the resolution plays a big role here. The smaller region
that a bead represents, the fewer contacts are found in Hi-C data. In
other words, unless fine grained simulations are run, it is likely
that some contacts have to be discarded.

The scaling targets the bead with most contacts, and is hence
sensitive to outliers. Possibly it might be better to scale to the
average number of contacts. If the Hi-C data is normalized so that
each bead has the same number of contacts, this comment does not apply.

Normalization with SKA gives each bead the same number of
contacts. That does to some extent hide any differences in
accessibility between regions. On the other hand it hides differences
in compactness. It is an open question how and if balancing should be
used in this context.

Please note that there are other, possibly better ways to convert Hi-C
matrices to contact probability matrices.

## REQUIRED ARGUMENTS
**\--hFile file.double**
: The Hi-C matrix, encoded as raw double.

**\--lFile file.uint8**
: An array with chromosome labels encoded as uint_8.

**\--nCont n**
: Set the max number of contacts per bead/bin.

**\--nStruct n**
: The number structures that will be generated. **hic2cpm** needs to
  know this since it affects which contacts will be used or not. For
  example a contact with probability 0.01 will be enabled only with
  100 or more structures.

**\--aOut file.double**
: Specify the output file name.


## OPTIONAL ARGUMENTS
**\--mode_eq**
: Enable balancing (KR normalization). This will force each bead to
  have the same number of contacts.

**\--y**
: Keep chrY encoded by 24. Without this argument chromosome Y will be removed.

**\--usage, \--help**
: Show a brief summary of available arguments.

# SPRITE2CPM
Can be used to convert sprite data (in .cluster files) to contact
probability matrices suitable for chromflock.

## Required arguments:
**\--sFile file**
: sprite file (.cluster)

## Optional arguments
**\--oFile file**
: output file

**\--binSize N**
: size of each bin

**\--balance**
: enable matrix balancing (see notes in hic2cpm)

**\--maxSize s**
: max cluster size, any cluster larger than that will be ignored.

**\--minSize s**
: min cluster size, any cluster smaller than that will be ignored

**\--printBins**
: Generate a .bins file where genomic coordinates are converted to bin number.

**\--hic**
: Parse FAHiC-files

**\--chr c**
: Only extract chr c

**\--verbose N**
: Set verbose level to N, default=0

**\--version**
: Show version

# STRING2ANY and ANY2STRING
The two commands can be used to quickly create or inspect data to/
written from chromflock.

## EXAMPLES

``` shell
$ chromflock string2any L.uint8 uint8_t 1 2 3
$ chromflock any2string uint8_t L.uint8
1
2
3
$ chromflock string2any C.double double 1.2 2.4
$ chromflock any2string double C.double
1.200000
2.400000
```


# SEE ALSO
**mflock**, **aflock**

# WEB PAGE
[https://github.com/elgw/chromflock/](https://github.com/elgw/chromflock/)

# REPORTING BUGS
Please report bugs at
[https://github.com/elgw/chromflock/issues/](https://github.com/elgw/chromflock/issues/)

# COPYRIGHT
Copyright © 2022 Erik Wernersson.  License GPLv3+: GNU GPL version 3
or later <https://gnu.org/licenses/gpl.html>.  This is free software:
you are free to change and redistribute it.  There is NO WARRANTY, to
the extent permitted by law.
