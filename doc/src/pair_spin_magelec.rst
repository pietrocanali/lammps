.. index:: pair_style spin/magelec

pair_style spin/magelec command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style spin/magelec cutoff

* cutoff = global cutoff pair (distance in metal units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style spin/magelec 4.5
   pair_coeff * * magelec 4.5 0.00109 1.0 1.0 1.0

The legacy coefficient form is

.. code-block:: LAMMPS

   pair_coeff I J magelec cutoff U ux uy uz

An instantaneous polarization can instead drive the interaction:

.. code-block:: LAMMPS

   pair_coeff I J polarization cutoff gamma center Ccore Cshell &
              neighbors Nspec core1 shell1 count1 weight1 ... polcut Rpol

For the Bi-centered BiFeO3 convention used in the accompanying development
test, Fe shells are magnetic type 4 and the coefficient is

.. code-block:: LAMMPS

   pair_coeff 4 4 polarization 4.5 1.09e-4 center 2 5 &
              neighbors 2 1 4 8 0.125 3 6 12 0.25 polcut 4.0

To couple each magnetic shell to its own core-centered local cell, use
*polarization/local* with the same argument structure.  For Fe-centered
cells in periodic rhombohedral BiFeO3:

.. code-block:: LAMMPS

   pair_coeff 4 4 polarization/local 4.5 1.09e-4 center 1 4 &
              neighbors 2 2 5 8 0.125 3 6 6 0.5 polcut 4.0

Here types 1/4, 2/5 and 3/6 are Fe, Bi and O core/shell pairs, respectively.
The Fe cell includes its own shell, eight Bi pairs with weight 1/8 and six
O pairs with weight 1/2: 29 contributing particles.  Changing the center
alone while keeping the Bi-centered oxygen counts and weights does not
define the Fe-centered cell.

Description
"""""""""""

Style *spin/magelec* computes a magneto-electric interaction between
pairs of magnetic spins. According to the derivation reported in
:ref:`(Katsura) <Katsura1>`, this interaction is defined as:

.. math::

   \vec{\omega}_i & = -\frac{1}{\hbar} \sum_{j}^{Neighb} \vec{s}_{j}\times\vec{D}(r_{ij}) \\
   \vec{F}_i & = -\sum_{j}^{Neighb} \frac{\partial D(r_{ij})}{\partial r_{ij}} \left(\vec{s}_{i}\times \vec{s}_{j} \right) \cdot \vec{r}_{ij}

where :math:`\vec{s}_i` and :math:`\vec{s}_j` are neighboring magnetic
spins of two particles.

From this magneto-electric interaction, each spin i will be submitted
to a magnetic torque omega, and its associated atom can be submitted to a
force F for spin-lattice calculations (see :doc:`fix nve/spin <fix_nve_spin>`),
such as:

.. math::

    \vec{F}^{i} & = -\sum_{j}^{Neighbor} \left( \vec{s}_{i}\times \vec{s}_{j} \right) \times \vec{E} \\
    \vec{\omega}^{i} = -\frac{1}{\hbar} \sum_{j}^{Neighbor} \vec{s}_j \times \left(\vec{E}\times r_{ij} \right)

with h the Planck constant (in metal units) and :math:`\vec{E}` an
electric polarization vector.  The norm and direction of E are giving
the intensity and the direction of a screened dielectric atomic
polarization (in eV).

The legacy *magelec* coefficient form retains this historical input-vector
normalization and its numerical behavior.

With either polarization coefficient form the electric vector is recomputed
at every pair-force evaluation and before each fixed-position spin-update
block of :doc:`fix nve/spin <fix_nve_spin>`, including the block after the
position update.  This corrects the stale polarization used by earlier
development versions during moving-lattice spin integration.  For each owned center core, the center shell
having the same molecule ID is selected, followed by the requested number of
nearest cores for every neighbor specification.  Each selected core is paired
with the shell having the same molecule ID.  Triclinic minimum-image
displacements are represented by LAMMPS's translated ghost coordinates; the
center's own periodic images are excluded, and only the nearest image of every
other atom ID is retained.  Every selected core and shell is included with
the specified per-particle weight.  The cell dipoles are reduced over MPI ranks:

.. math::

   \mathbf P = \frac{16.02176634}{V}
     \sum_{c}\sum_{a\in c} w_a q_a(\mathbf r_a-\mathbf r_c)

Here charges are in units of the electron charge, distances are in Angstrom,
the instantaneous box volume :math:`V` is in Angstrom cubed, and the resulting
polarization is in C/m\ :sup:`2`.  The conversion uses
:math:`1\ e/\mathrm{Angstrom}^2=16.02176634\ \mathrm{C/m^2}`.  The dynamic
Hamiltonian is

.. math::

   H_\mathrm{ME}=\gamma\sum_{i<j}
     (\mathbf P\times\mathbf e_{ij})\mathbin{\cdot}(\mathbf s_i\times\mathbf s_j),

where :math:`\mathbf e_{ij}` points from atom *i* to atom *j*.  *gamma* is in
eV per C/m\ :sup:`2`; unlike the legacy input vector, :math:`\mathbf P` is not
normalized.

For *polarization/local*, each cell instead supplies

.. math::

   \mathbf P_i &= \frac{16.02176634 N_c}{V}
      \sum_{a\in c(i)}w_aq_a(\mathbf r_a-\mathbf r_{c(i)}), \\
   \mathbf P_{ij} &= \frac{\mathbf P_i+\mathbf P_j}{2}, \\
   H_\mathrm{ME} &= \sum_{i<j}\gamma_{ij}
      (\mathbf P_{ij}\times\mathbf e_{ij})\mathbin{\cdot}(\mathbf s_i\times\mathbf s_j).

Here :math:`c(i)` is the core sharing the magnetic shell's molecule ID and
:math:`N_c` is the global number of center cores.  The volume assigned to
each cell is the instantaneous mean :math:`V/N_c`, not a geometric local
volume.  Polarization magnitudes and directions fluctuate freely.  The
bond direction is still measured between the spin-bearing shells.

The local gradient is

.. math::

   \mathbf G_i = \frac{\partial H}{\partial\mathbf P_i}
     = \frac12\sum_j\gamma_{ij}\mathbf e_{ij}\times(\mathbf s_i\times\mathbf s_j).

Every contributor receives :math:`-16.02176634(N_c/V)w_aq_a\mathbf G_i`,
and the center receives the opposite force, in addition to the separate
bond-direction forces.  Cell polarizations and gradients are communicated
between core owners and spin-shell owners; these pairs need not be on the
same MPI rank.  The full spin list supplies half-weighted directed energies
and bond forces, while the spin derivative is not halved.

The dynamic form differentiates this same Hamiltonian with respect to the
spins, the pair direction, every polarization contributor, and the explicit
:math:`1/V` factor.  It therefore supplies magnetic torques, reciprocal
atomic forces when ``fix nve/spin lattice moving`` is used, and the associated
six-component virial (including the diagonal volume-normalization term).
Membership is reselected from current coordinates at each polarization
evaluation using the available neighbor list, so periodic image crossings
and triclinic boxes are handled without storing wrapped reference positions.

The current polarization is exposed as the three-element pair extra vector:

.. code-block:: LAMMPS

   compute pol all pair spin/magelec
   thermo_style custom step c_pol[1] c_pol[2] c_pol[3]

In local mode this remains the mean :math:`N_c^{-1}\sum_i\mathbf P_i`.
The scalar compute is the interaction energy; it cannot in general be
reconstructed using only the mean vector.  The explicit volume-normalization
virial is included globally and partitioned equally between the spin-bond
endpoints for per-atom stress.

For a hybrid pair style the same compute syntax selects the
``spin/magelec`` sub-style.  Dynamic cell-definition data and pair
coefficients are written to restart files; polarization itself is transient
and is recomputed during the first force evaluation after restart.
The local mode uses an additional value in the existing restart coefficient
field, so this implementation can still read the preceding development
version's fixed and global-polarization restart data.

More details about the derivation of these torques/forces are reported in
:ref:`(Tranchida) <Tranchida4>`.

----------

Restrictions
""""""""""""

All the *pair/spin* styles are part of the SPIN package.  These styles
are only enabled if LAMMPS was built with this package, and if the
atom_style "spin" was declared.  See the :doc:`Build package <Build_package>` page for more info.

The dynamic *polarization* form additionally requires per-atom charges,
positive atom and molecule IDs, and explicit core/shell type pairs.  Thus an
atom style such as ``atom_style hybrid full spin`` is appropriate for a
core-shell model.  Every center must have one center shell and the requested
number of core-shell neighbors within *Rpol*.  Missing membership is a fatal
error that reports the center atom ID, expected and observed counts, and the
polarization cutoff.  All type numbers and both cutoffs must be valid; one
cell definition is shared by all dynamic coefficients of this pair-style
instance.

The local mode additionally requires full periodicity, neutral weighted
cells (absolute charge tolerance :math:`10^{-8}` electron charges), distinct
core/shell species types, positive finite weights and unique core-shell
partners.  Both pair-coefficient types must equal the center-shell type,
and those shells must have positive spin magnitude.  A dynamic pair cutoff
must not exceed the pair-style cutoff.  The communicated ghost range is
automatically enlarged when the polarization cutoff is larger than the
magnetic cutoff.

This local prescription is intended for bulk perovskite structures with
intact coordination, not surfaces, vacuum regions, vacancies or reconstructed
cells.  Supercells must contain enough distinct atom IDs to supply every
cell neighbour: repeated periodic images of one atom cannot fill several
cell sites.  Missing membership reports the center, type, count and cutoff.
The polarization cutoff must include both members of every selected pair,
with a margin for thermal displacements and strain.

Nearest-neighbour membership and magnetic cutoffs are discrete.  The force
derivatives apply between membership changes; exchanging neighbours or
crossing a sharp magnetic cutoff can cause energy discontinuities.  Check
coordination and cutoff margins for the intended trajectories.

Fe-centered and Bi-centered local polarizations are different microscopic
definitions even in the R3c phase.  A [111] polar direction does not make a
Bi-to-Fe shift a crystal translation.  The scalar *gamma* remains a
phenomenological parameter; this term alone does not add a separate
octahedral-rotation-mediated magnetic coupling or magnetic anisotropy.

Related commands
""""""""""""""""

:doc:`atom_style spin <atom_style>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`pair_style spin/exchange <pair_spin_exchange>`, :doc:`pair_eam <pair_eam>`,

Default
"""""""

none

----------

.. _Katsura1:

**(Katsura)** H. Katsura, N. Nagaosa, A.V. Balatsky. Phys. Rev. Lett., 95(5), 057205. (2005)

.. _Tranchida4:

**(Tranchida)** Tranchida, Plimpton, Thibaudeau, and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).
