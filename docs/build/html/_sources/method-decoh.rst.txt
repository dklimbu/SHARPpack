Decoherence
===========

Energy-based decoherence correction :cite:`decoh2007` applies exponential damping to inactive adiabatic states, while renormalizing the active state to preserve the total wavefunction norm (:math:`\sum |c_i|^2 = 1`).

.. math::

   c_i \rightarrow c_i \exp\left(-\frac{\Delta t}{\tau_{ia}}\right),  \quad \forall \; i \ne a,

.. math::

   c_a = c_a \left[\frac{1 - \sum_{i \ne a}|c_i|^2}{|c_a|^2}\right]^{1/2},

with

.. math::

   \tau_{ia} = \frac{\hbar}{|E_i - E_a|}
   \left( C + \frac{E_0}{T_a} \right),

where

- :math:`E_i` is the energy of inactive state :math:`i`,
- :math:`E_a` is the energy of active state :math:`a`,
- :math:`T_a` is the kinetic energy on the active surface,
- :math:`C` and :math:`E_0` are empirical constants.


