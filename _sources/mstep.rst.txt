Molecular dynamics (M-step)
---------------------------

The m-step uses simulated annealing with a decreasing temperature.
Currently 8 temperatures are used: :math:`T=.7^s`, where :math:`s =
0, ..., 6` and then finally :math:`T=0`.

If chromosome compression is enabled, it is turned off when :math:`T=0`
in order for the structure to relax.

The forces in common with Alber are:

-  :math:`F_v`, volume exclusion, i.e. beads are not allowed to overlap.

-  :math:`F_s`, confinement in the sphere.

-  :math:`F_a`, attraction between beads in contacts.

-  :math:`F_c`, chromosomal compression.

As well as:

-  :math:`F_r`, radial positioning (GPSeq).

-  :math:`F_b`, Brownian force.

-  :math:`F_v`, viscosity or dampening proportional to the velocity of
   each bead.

Definition of forces
~~~~~~~~~~~~~~~~~~~~

Here follows a list of the forces that chromflock uses. The Brownian
force:

.. math:: F_b(i) = c_bs

where :math:`s` is drawn from an isotropic 3D Gaussian with
:math:`\sigma=1`. The drag force (viscosity):

.. math:: F_v(i) = -\eta v(i)

The other forces are defined in terms of their errors (or potential). We
set :math:`d_{ij} = ||X_i - X_j||`, and :math:`r_i=||X_i||`, then the
volume exclusions potential that keeps the beads from overlapping is:

.. math::

   E_v(i,j) = \begin{cases}
       c_v(d_{ij} - 2R_b)^2 & \text{if } d_{ij} < 2R_b,\\
       0 & \text{if } d_{ij} > 2R_b.
     \end{cases}

The potential that keeps the beads inside the domain (nuclei):

.. math::

   E_d(i) = \begin{cases}
       c_d(r_i+R_b-R_s)^2 & \text{if } r_i > R_s-R_b, \\
       0 & \text{if } r_i < 1-R_b
     \end{cases}
     \label{eq:domError}

And the potential that keeps beads attracted to each other:

.. math::

   E_a(i,j) = \begin{cases}
       c_a(d_{ij}-R_c)^2 & \text{if } W_{ij} =1 \mbox{ and } d_{ij}>R_c, \\
       0 & \text{if } W_{ij} =0.
     \end{cases}

And the potential for radial preference:

.. math::

   E_r(i) = \begin{cases}
       c_r (r_i-g_i)^2 & \text{if } g_i \text{is finite} \\
       0 & \text{else}
     \end{cases}

.. math:: E_c(i) = c_c(p_i-m_k)^2

when bead :math:`i` belongs to chromosome :math:`k`. We let the volume
exclusion force vary with time as:

.. math::

   F_v(x) = \frac{1}{2}\left(1+\mbox{erf}(\beta(p-.5)\right)
     \label{eq:kvol}

where we let :math:`\beta`\ =5.

The total error is

.. math:: E = \sum_{i=1}^N \left(E_d(i) + E_c(i) + E_r(i)\right) + \sum_{i,j} \left(E_v(i,j)+E_a(i,j)\right)

and hence the total forces are:

.. math:: F = \nabla E + F_b + F_v

It is inefficient to calculate :math:`\nabla E` by finite differences,
and hence we have derived analytic expressions which we have validated
against the numerical gradient.

Time evolution
~~~~~~~~~~~~~~

We use Verlet integration [1]_ to
evolve the beads in time:

.. math:: X_{t+1} = 2X_t - X_{t-1} - (\Delta t)^2F

Parameters
~~~~~~~~~~

Default parameters are: :math:`R_s = 1`, :math:`V_q = 0.2`
:math:`R_c = (2.1+0.9p)R_b` :math:`\eta = 0.5`, :math:`\Delta t = 0.15`,
:math:`c_a = 1`, :math:`c_b=1`. :math:`c_c=1`, :math:`c_d=1`,
:math:`c_r=0` or :math:`0.001` if the radial force is enabled,
:math:`c_s=1`, :math:`c_v=1`. Where :math:`p\in[0,1]` is the proportion
of steps taken.

Finding neighbours
~~~~~~~~~~~~~~~~~~

An indirect hash table constructed with count-sort is used. It has low
memory usage since the load factor is 1. It might be faster to use open
addressing with linear probing since there are limits for how many beads
there can be per bucket due to self-repulsion.

References
~~~~~~~~~~

.. [1] T. Schlick (2010) Molecular Modeling and Simulation: An Interdisciplinary Guide, Spinger `doi <https://doi.org/10.1007/978-1-4419-6351-2>`_
