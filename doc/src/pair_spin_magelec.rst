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

With the *polarization* coefficient form the electric vector is recomputed
at every pair-force evaluation.  For each owned center core, the center shell
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

The dynamic form differentiates this same Hamiltonian with respect to the
spins, the pair direction, every polarization contributor, and the explicit
:math:`1/V` factor.  It therefore supplies magnetic torques, reciprocal
atomic forces when ``fix nve/spin lattice moving`` is used, and the associated
six-component virial (including the diagonal volume-normalization term).
Membership is rebuilt with the neighbor lists, so periodic image crossings
and triclinic boxes are handled without storing wrapped reference positions.

The current polarization is exposed as the three-element pair extra vector:

.. code-block:: LAMMPS

   compute pol all pair spin/magelec
   thermo_style custom step c_pol[1] c_pol[2] c_pol[3]

For a hybrid pair style the same compute syntax selects the
``spin/magelec`` sub-style.  Dynamic cell-definition data and pair
coefficients are written to restart files; polarization itself is transient
and is recomputed during the first force evaluation after restart.

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
