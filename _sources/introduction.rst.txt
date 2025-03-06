Chromflock
==========

chromflock is a software package to deconvolve bulk Hi-C data into
putative single-cell structures. The overall scheme is borrowed from
PGS [1]_ from Frank Alber lab although it is not a direct fork.

Chromflock will:

- Work with haploid as well as diploid structures.

- Can integrate GPSeq data for radial preferences.

- Supports spherical as well as an ellipsoidal domain for the beads.

The documentation is not complete at the moment and it is suggested
that anyone interested in chromflock start by reading the papers from
Frank Alber lab, since much of the terminology used here is introduced
by them.

Input and output
----------------

As input, chromflock takes:

- A contact probability matrix, :math:`A`, where :math:`A_{ij}`
  denotes the proportion of structures that should have a contact
  between site :math:`i` and :math:`j`, i.e. a number in
  :math:`[0,1]`.

-  A label matrix, :math:`L`, which points out which chromosome each
   bead belongs to.

-  A vector specifying radial preferences for each bead, :math:`R`.
   Possibly together with a vector specifying the probability of each
   radial preference to be used.

-  The number of structures to generate, :math:`N`.

-  A list of theta values, :math:`\theta_1 (=1) > \theta_2 > ... > \theta_T`,
   where :math:`T-1` is the number of epochs to be run and :math:`\theta_t`
   specifies which contacts from the :math:`A`-matrix that
   should be used in epoch :math:`t`.

-  The number of times that each epoch should be re-run, :math:`C`.

-  A lua script which defines the dynamics of the simulations.

The output consists of:

-  Coordinates for all beads in all structures, :math:`X`.

-  The contact indicator tensor, :math:`W^{(n)}_{i,j}` which says which
   contacts are assigned to which structures.

-  Log files and some pre-calculated summary statistics such
   as the average radial profile.

Workflow and time complexity
----------------------------


-  Initialize :math:`S` structures with :math:`N` beads. (aflock)

-  Run the M-step to find initial structures. (mflock)

-  For each of the :math:`(T-1)` epochs:

   -  For each of the :math:`C` AM-cycles:

      -  (re)-assign contacts (aflock)

      -  Find structures (mflock)

-  Finalize (aflock)

The time consumption is typically dominated by the M-step (mflock) and
the number of mflock commands to be run is :math:`N_M = S((T-1)C+1)`.

Typical
values would be :math:`S=10,000`, :math:`T=6`, :math:`C=4`, i.e.
:math:`N_M=210,000`. Using about :math:`3,000` beads, an M-step takes
:math:`12` s. Using a machine with 40 cores this will take:
:math:`210,000\times12/(40\times3,600)` s or about :math:`= 17.5` hours.

Geometry
--------

The domain of the simulations is always an ellipsoid, however for
simplicity we present it for a spherical domain first. In that case the
sphere has a radius :math:`R_s=1`. The beads have a radius, :math:`R_b`
and beads in contact are set to have a distance below
:math:`R_c=4.0R_b`.

Typically :math:`R_b` is set to get a certain volume occupancy. The
volume of the spherical simulation domain is:

.. math:: V_s = \frac{4\pi}{3}

and the volume of the :math:`N` beads is

.. math:: V_b = N(R_b)^3\frac{4\pi}{3}

setting :math:`V_b/V_s=V_q` gives

.. math:: R_b = \sqrt[3]{\frac{V_q}{N}}.

See TODO for the elliptical geometry.


.. [1] Hua, N., Tjong, H., Shin, H. et al. Producing genome structure
   populations with the dynamic and automated PGS software. Nat Protoc
   13, 915–926 (2018). https://doi.org/10.1038/nprot.2018.008
   https://github.com/alberlab/pgs
