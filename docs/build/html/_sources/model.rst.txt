Models
======

.. _tab-model:

.. list-table:: Models
   :widths: 20 80
   :header-rows: 1

   * - **Keyword**
     - **Description**
   * - `tully1`
     - Tully's Simple Avoided Crossing Model
   * - `tully2`
     - Tully's Dual Avoided Crossing Model
   * - `tully3`
     - Tully's Extended Coupling with Reflection Model
   * - `morse1`
     - Morse Model I – photo-dissociation
   * - `morse2`
     - Morse Model II – photo-dissociation
   * - `morse3`
     - Morse Model III – photo-dissociation
   * - `db2lchain`
     - 2-state model coupled with linear chain
   * - `db3lchain`
     - 3-state super-exchange model with linear chain coupling
   * - `spinboson`
     - Spin-Boson model coupled to N-harmonic bath

Tully Models
------------

Diabatic Hamiltonians:

**Model I**

.. math::

   V_{11}(R) = 
   \begin{cases}
      A(1 - e^{-B R}) & \text{if } R > 0 \\
     -A(1 - e^{B R})  & \text{if } R < 0
   \end{cases}

.. math::
   V_{22}(R) = -V_{11}(R)

.. math::
   V_{12}(R) = V_{21}(R) = C e^{-D R^2}

**Model II**

.. math::
   V_{11}(R) = 0

.. math::
   V_{22}(R) = -A e^{-B R^2} + E_0

.. math::
   V_{12}(R) = V_{21}(R) = C e^{-D R^2}

**Model III**

.. math::
   V_{11}(R) = A, \quad V_{22}(R) = -A

.. math::

   V_{12}(R) = V_{21}(R) = 
   \begin{cases}
     B e^{C R}, & \text{if } R < 0 \\
     B(2 - e^{-C R}), & \text{if } R > 0
   \end{cases}

**Parameters (atomic units):**

.. list-table:: Tully Model Parameters
   :widths: 10 15 15 15
   :header-rows: 1

   * - **Parameter**
     - **Model I**
     - **Model II**
     - **Model III**
   * - A
     - 0.01
     - 0.1
     - 6e-4
   * - B
     - 1.6
     - 0.28
     - 0.1
   * - C
     - 0.005
     - 0.015
     - 0.9
   * - D
     - 1.0
     - 0.06
     - --
   * - E₀
     - --
     - 0.05
     - --

Morse Potential Models
----------------------

.. math::

   V_{ii}(R) = D_{e,i} \left(1 - e^{-\beta_i(R - R_{e,i})} \right)^2 + c_i

.. math::

   V_{ij}(R) = A_{ij} e^{-a_{ij}(R - R_{ij})^2}, \quad V_{ji}(R) = V_{ij}(R)

.. list-table:: Numerical values of the parameters related to the three Morse potential models (in a.u.)
   :widths: 10 10 10 10 10 10 10 10 10 10
   :header-rows: 2
   :align: center

   * - 
     - 
     - **Model I**
     - 
     - 
     - **Model II**
     - 
     - 
     - **Model III**
     - 
   * - **Param**
     - 1
     - 2
     - 3
     - 1
     - 2
     - 3
     - 1
     - 2
     - 3
   * - :math:`D_{e,i}`
     - 0.003
     - 0.004
     - 0.003
     - 0.02
     - 0.01
     - 0.003
     - 0.02
     - 0.02
     - 0.003
   * - :math:`\beta_i`
     - 0.65
     - 0.60
     - 0.65
     - 0.65
     - 0.40
     - 0.65
     - 0.40
     - 0.65
     - 0.65
   * - :math:`R_{e,i}`
     - 5.0
     - 4.0
     - 6.0
     - 4.5
     - 4.0
     - 4.4
     - 4.0
     - 4.5
     - 6.0
   * - :math:`c_i`
     - 0.00
     - 0.01
     - 0.006
     - 0.00
     - 0.01
     - 0.02
     - 0.02
     - 0.00
     - 0.02
   * - :math:`ij`
     - 12
     - 13
     - 23
     - 12
     - 13
     - 23
     - 12
     - 13
     - 23
   * - :math:`A_{ij}`
     - 0.002
     - 0.0
     - 0.002
     - 0.005
     - 0.005
     - 0.0
     - 0.005
     - 0.005
     - 0.0
   * - :math:`R_{ij}`
     - 3.40
     - 0.0
     - 4.80
     - 3.66
     - 3.34
     - 0.0
     - 3.40
     - 4.97
     - 0.0
   * - :math:`a_{ij}`
     - 16.0
     - 0.0
     - 16.0
     - 32.0
     - 32.0
     - 0.0
     - 32.0
     - 32.0
     - 0.0

Linear Chain Models
----------------------


Spin–Boson Model
----------------

The Spin–Boson Hamiltonian :cite:`Leggett1987` is

.. math::

   H = H_s + H_b + H_{sb}.

Two–level system Hamiltonian:

.. math::

   H_s = \epsilon \sigma_z + \Delta \sigma_x =
   \begin{pmatrix}
   \epsilon & \Delta \\
   \Delta & -\epsilon
   \end{pmatrix},

with bias :math:`\epsilon` and coupling :math:`\Delta`.

Bath Hamiltonian:

.. math::

   H_b = \sum_j \left( \frac{P_j^2}{2 M_j}
   + \frac{1}{2} M_j \omega_j^2 R_j^2 \right).

System–bath coupling:

.. math::

   H_{sb} = \sigma_z \sum_j c_j R_j =
   \begin{pmatrix}
   \sum_j c_j R_j & 0 \\
   0 & -\sum_j c_j R_j
   \end{pmatrix}.

----

Spectral Density
----------------

General definition for :math:`N` bath modes :cite:`Miller1999-debye, Templaar2018`:

.. math::

   J(\omega) = \frac{\pi}{2} \sum_{j=1}^{N}
   \frac{c_j^2}{M_j \omega_j} \, \delta(\omega - \omega_j),

where,  
  - :math:`M_j` = mass of :math:`j`-th oscillator,  
  - :math:`c_j` = coupling constant,  
  - :math:`\omega_j` = oscillator frequency.

Debye form
^^^^^^^^^^

.. math::

   J_D(\omega) = \frac{E_r}{2} \,
   \frac{\omega \, \omega_c}{\omega^2 + \omega_c^2},

where,
  - :math:`E_r` = reorganization energy,  
  - :math:`\omega_c` = cutoff frequency.

----

**Debye Bath Discretization**

Frequencies and couplings are discretizeda as :cite:`Wang:2001_bath_discret,Hanna:2013_bath_discret`:

.. math::

   \omega_j = \omega_c \cdot
   \tan\left( \frac{j}{N} \tan^{-1}\left( \frac{\omega_{\text{max}}}{\omega_c} \right) \right),

.. math::

   c_j = \omega_j \sqrt{
   \frac{M_j E_r}{\pi N} \,
   \tan^{-1}\left( \frac{\omega_{\text{max}}}{\omega_c} \right)
   },

with :math:`\omega_{\text{max}} = 20 \omega_c`.

----

Ohmic form
^^^^^^^^^^

.. math::

   J_O(\omega) = \frac{\pi}{2} \, \alpha \, \omega \, e^{\omega/\omega_c},

where,
  - :math:`\alpha` = Kondo parameter that controls the strength of the coupling,  
  - :math:`\omega_c` = cutoff frequency.

**Ohmic Bath Discretization**

Frequencies and couplings are discretizeda as:

.. math::

   \omega_j = \omega_c \cdot
   \log\left(1 - \frac{j}{1+N} \right),

.. math::

   c_j = \omega_j \sqrt{
   \frac{M_j \alpha}{1+N}},

