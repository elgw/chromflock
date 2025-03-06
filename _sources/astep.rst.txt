Contact assignment (A-step)
---------------------------

The A-step runs in several modes:

Initialization:
^^^^^^^^^^^^^^^

Each structure, :math:`n` is assigned only the contacts where :math:`A=1`, i.e.

.. math::

   W_{ij}^{(n)} = \begin{cases} 1 & \text{if } A_{ij} = 1 \\ 0 &
       \text{if } A_{ij} < 1 \\
     \end{cases}

Update:
^^^^^^^

Using an interval :math:`\theta_l \geq \theta \geq \theta_h`.

-  Load :math:`W`, i.e., :math:`W^{(n)}`, for all structures
   :math:`n\in \{1, ..., N\}`.

-  For each pair of beads, :math:`(i,j)`, :math:`i\neq j`: Let :math:`a`
   = :math:`A_{i,j}`. If :math:`a \in [\theta_l, \theta_h[` then
   :math:`s = \lfloor a \times N \rfloor` structures should have a contact between
   :math:`i` and :math:`j`. Lets define
   :math:`d=||X^{(n)}_i - X^{(n)}_j||`. Then

   .. math::

      D_n = \begin{cases}
              d_n + 2 & \text{if } W^{(n)}_{i,j} = 1 \text{ and } d_n > d_c \\
              d_n & \text{otherwhise}. \\
            \end{cases}

-  :math:`W^{(n)}_{i,j}` is set to :math:`0` for all :math:`n`. Then the
   :math:`s` structures for which :math:`D_n` is smallest are assigned a
   contact between :math:`i` and :math:`j`, i.e., :math:`W_{i,j}` is set
   to :math:`1` for them.

Final:
^^^^^^

The average radial profile, the contact map and the number of contact
violations are calculated.
