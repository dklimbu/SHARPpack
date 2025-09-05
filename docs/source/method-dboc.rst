Diagonal Born–Oppenheimer Correction (DBOC)
===========================================

When the nuclear kinetic energy operator

.. math::

   \hat{T}_\mathrm{n} = -\sum_{I} \frac{\hbar^2}{2 M_I} \nabla_{\mathbf{R}_I}^2

acts on the adiabatic electronic wavefunctions both diagonal and off-diagonal non-adiabatic couplings arise. :cite:`dboc2016` The diagonal term behaves like a potential energy correction and modifies the adiabatic
PESs, :math:`V_\alpha(\mathbf{R})`, as

.. math::

   \tilde{V}_{\alpha}(\mathbf{R}) =
      V_\alpha(\mathbf{R}) + V_{\mathrm{DBOC}}^{(\alpha)}(\mathbf{R}).

Here, :math:`\mathbf{R} \equiv \{\mathbf{R}_1, \mathbf{R}_2, \dots, \mathbf{R}_n\}`
denotes the ring-polymer coordinates containing *n* beads.
The correction term is

.. math::

   V_{\mathrm{DBOC}}^{(\alpha)}(\mathbf{R}) =
      \sum_{I=1}^N \sum_{k=1}^n
      \frac{\hbar^2}{2M_{I,k}}
      \sum_{\gamma \neq \alpha}
      \left| \mathbf{d}_{I,k,\alpha\gamma}(\mathbf{R}_k) \right|^2.

----

Total force on bead *k* for adiabatic state :math:`\alpha` is

.. math::

   \mathbf{F}_{I,k}
      = -\nabla \tilde{V_\alpha}(\mathbf{R}_k)
      = -\nabla V_\alpha(\mathbf{R}_k)
        - \nabla V_{\mathrm{DBOC}}^{(\alpha)}(\mathbf{R}_k)
      = \mathbf{F}_{I,k}^{(0)} + \Delta \mathbf{F}_{I,k},

where the DBOC force contribution is

.. math::

   \Delta \mathbf{F}_{I,k} =
      \sum_{\alpha \ne \gamma}
      \frac{1}{M_{I,k}}
      \mathbf{d}_{I,k,\alpha \gamma}(\mathbf{R}_k)
      \cdot
      \nabla \mathbf{d}_{I,k,\alpha \gamma}(\mathbf{R}_k).

The non-adiabatic coupling vectors (NACVs) are evaluated at each bead *k*:

.. math::

   \mathbf{d}_{I,k,\alpha \gamma} =
      \frac{\langle\psi_{k,\alpha}|
         \nabla_{R_{I,k}} H(\mathbf{R}_k)|\psi_{k,\gamma}\rangle}
           {E_\gamma(\mathbf{R}_k) - E_\alpha(\mathbf{R}_k)}.

----

DBOC-corrected forces are computed numerically via finite differences:

.. math::

   \nabla \mathbf{d}_{I,k,\alpha \gamma} =
      \frac{
        \mathbf{d}_{I,k,\alpha \gamma}(\mathbf{R}_k+\delta)
        - \mathbf{d}_{I,k,\alpha \gamma}(\mathbf{R}_k-\delta)
      }{2\delta}.

.. warning::
   The displacement parameter :math:`\delta` must be chosen to ensure numerical stability and convergence of the DBOC forces during dynamics simulations.

