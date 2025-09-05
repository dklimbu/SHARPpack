Models
======

There are several model systems implemented in **SHARP Pack**, which are summarized in the following table.

.. _tab-model:

.. list-table:: Models
   :widths: 20 80
   :header-rows: 1

   * - **Keyword**
     - **Description**
   * - `tully1`
     - Tully Model-I: Simple avoided crossing model :cite:`Tully:1990`
   * - `tully2`
     - Tully Model-II: Dual avoided crossing model :cite:`Tully:1990`
   * - `tully3`
     - Tully Model-III: Extended coupling with reflection model :cite:`Tully:1990`
   * - `morse1`
     - Morse Model-I: photo-dissociation model :cite:`morse_model`
   * - `morse2`
     - Morse Model-II: photo-dissociation model :cite:`morse_model`
   * - `morse3`
     - Morse Model-III: photo-dissociation model :cite:`morse_model`
   * - `db2lchain`
     - 2-state model coupled with linear chain model :cite:`Limbu:2023,Parandekar:2005`
   * - `db3lchain`
     - 3-state super-exchange model with linear chain model :cite:`Sifain:2016`
   * - `superexchange`
     - 3-state super-exchange model :cite:`superexchange2014`
   * - `spinboson`
     - Spin-Boson model coupled to N-harmonic bath :cite:`Leggett1987`
   * - `dboc1`
     - Flat Born-Oppenheimer PES model :cite:`pssh2009,dboc2016`

Tully Models
------------

Diabatic Hamiltonians and parameters for Tully models :cite:`Tully:1990` are:

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

.. list-table:: Tully Model Parameters (in a.u.)
   :widths: 10 15 15 15
   :header-rows: 1

   * - **Parameter**
     - **Model I**
     - **Model II**
     - **Model III**
   * - A
     - 0.01
     - 0.1
     - 6 :math:`\times 10^{-4}`
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

Diabatic Hamiltonians and parameters for three Morse models :cite:`morse_model` are:

.. math::

   V_{ii}(R) = De_i \left(1 - e^{-\beta_i(R - Re_i)} \right)^2 + c_i

.. math::

   V_{ij}(R) = A_{ij} e^{-a_{ij}(R - R_{ij})^2}, \quad V_{ji}(R) = V_{ij}(R)


.. list-table:: Numerical values of the parameters related to the three Morse potential models (in a.u.)
   :widths: 15 12 12 12 12 12 12 12
   :header-rows: 1
   :align: center

   * - **Model-I**
     - :math:`De_i` 
     - :math:`\beta_i`
     - :math:`Re_i`
     - :math:`c_{i}`
     - :math:`A_{ij}`
     - :math:`R_{ij}`
     - :math:`a_{ij}`
   * - :math:`V_{11}`
     - 0.003
     - 0.65 
     - 5.0
     - 0.00
     - 
     - 
     - 
   * - :math:`V_{22}`
     - 0.004
     - 0.60 
     - 4.0
     - 0.01
     - 
     - 
     - 
   * - :math:`V_{33}`
     - 0.003
     - 0.65 
     - 6.0
     - 0.006
     - 
     - 
     - 
   * - :math:`V_{12}`
     - 
     -  
     - 
     - 
     - 0.002
     - 3.40
     - 16.0
   * - :math:`V_{23}`
     - 
     -  
     - 
     - 
     - 0.002
     - 4.80
     - 16.0
   * - **Model-II**
     - :math:`De_i`
     - :math:`\beta_i`
     - :math:`Re_i`
     - :math:`c_{i}`
     - :math:`A_{ij}`
     - :math:`R_{ij}`
     - :math:`a_{ij}`
   * - :math:`V_{11}`
     - 0.02
     - 0.65 
     - 4.5
     - 0.00
     - 
     - 
     - 
   * - :math:`V_{22}`
     - 0.01
     - 0.40 
     - 4.0
     - 0.01
     - 
     - 
     - 
   * - :math:`V_{33}`
     - 0.003
     - 0.65 
     - 4.4
     - 0.02
     - 
     - 
     - 
   * - :math:`V_{12}`
     - 
     -  
     - 
     - 
     - 0.005
     - 3.66
     - 32.0
   * - :math:`V_{13}`
     - 
     -  
     - 
     - 
     - 0.005
     - 3.34
     - 32.0
   * - **Model-III**
     - :math:`De_i`
     - :math:`\beta_i`
     - :math:`Re_i`
     - :math:`c_{i}`
     - :math:`A_{ij}`
     - :math:`R_{ij}`
     - :math:`a_{ij}`
   * - :math:`V_{11}`
     - 0.02
     - 0.40 
     - 4.0
     - 0.02
     - 
     - 
     - 
   * - :math:`V_{22}`
     - 0.02
     - 0.65 
     - 4.5
     - 0.00
     - 
     - 
     - 
   * - :math:`V_{33}`
     - 0.003
     - 0.65 
     - 6.0
     - 0.02
     - 
     - 
     - 
   * - :math:`V_{12}`
     - 
     -  
     - 
     - 
     - 0.005
     - 3.40
     - 32.0
   * - :math:`V_{13}`
     - 
     -  
     - 
     - 
     - 0.005
     - 4.97
     - 32.0


Linear Chain Model
------------------

Linear Chain Model is **two-level/three-level** system copuled to a signle atom in N-atom
linear chain, with the following anharmonic, nearest-neighbor potential energy function:

.. math::

   V(\mathbf{R}) = \sum_{i=1}^{N} V_M(R_i - R_{i+1})

where

.. math::

   V_M(R) = V_0 \left(a^2 R^2 - a^3 R^3 + 0.58 a^4 R^4 \right)

and :math:`R_{N+1}` is a fixed position and the atom farthest from the quantum systerm (:math:`R_N`) is connected to Langevin dynamics with friction constant :math:`\gamma`. :cite:`Limbu:2023,Parandekar:2005,Sifain:2016`

.. list-table:: Simulation Parameters for **two-level** Linear Chain Model 
   :widths: 20 20 80
   :header-rows: 1

   * - **Parameter**
     - **Value**
     - **Unit**
   * - :math:`N` 
     - 20
     - 
   * - :math:`m` 
     - 12.0
     - amu
   * - :math:`V_0` 
     - 175.0 
     - kJ/mol
   * - :math:`a` 
     - 4.0
     - :math:`\angstrom ^{-1}`
   * - :math:`\gamma` 
     - :math:`\text{10}^{14}`
     - :math:`s^{-1}`
   * - :math:`\Delta=\epsilon_2-\epsilon_1` 
     - 8.0
     - kJ/mol
   * - :math:`d_{12}` 
     - -6.0
     - :math:`\angstrom ^{-1}`

.. list-table:: Simulation Parameters for **three-level** Superexchange Linear Chain Model 
   :widths: 20 20 80
   :header-rows: 1

   * - **Parameter**
     - **Value**
     - **Unit**
   * - :math:`N` 
     - 20
     - 
   * - :math:`m` 
     - 12.0
     - amu
   * - :math:`V_0` 
     - 175.0 
     - kJ/mol
   * - :math:`a` 
     - 4.0
     - :math:`\angstrom ^{-1}`
   * - :math:`\gamma` 
     - :math:`\text{10}^{14}`
     - :math:`s^{-1}`
   * - :math:`\epsilon_1` 
     - 0
     - kJ/mol
   * - :math:`\epsilon_2` 
     - 39.0
     - kJ/mol
   * - :math:`\epsilon_3` 
     - 13.0
     - kJ/mol
   * - :math:`d_{12}` 
     - -6.0
     - :math:`\angstrom ^{-1}`
   * - :math:`d_{23}` 
     - 8.0
     - :math:`\angstrom ^{-1}`
   * - :math:`d_{13}` 
     - 0
     - :math:`\angstrom ^{-1}`


Super Exchange Model
--------------------

Diabatic Hamiltonians and parameters for 3-state super exchange model :cite:`superexchange2014` are defined as:

.. math::

   V_{ii}(R) = A_i

.. math::

   V_{ij}(R) = V_{ji}(R) = B_{ij}\, e^{-R^2/2}


.. list-table:: Super Exchange model parameter (in a.u.)
   :widths: 15 12 12 12 12 12 12
   :header-rows: 1
   :align: center

   * - 
     - :math:`i=1` 
     - :math:`i=2`
     - :math:`i=3`
     - :math:`ij=12`
     - :math:`ij=23`
     - :math:`ij=13`
   * - :math:`A_i`
     - 0
     - 0.01
     - 0.005
     -
     -
     -
   * - :math:`B_{ij}`
     - 
     -  
     - 
     - 0.001
     - 0.01
     - 0

Spin–Boson Model
----------------

The **spin-boson model** provides a theoretical framework for studying non-adiabatic dynamics in a condensed-phase environment. The corresponding Hamiltonian :cite:`Leggett1987` is given by:

.. math::

   H = H_s + H_b + H_{sb}.

Two–level system Hamiltonian:

.. math::

   H_s = \epsilon \sigma_z + \Delta \sigma_x =
   \begin{pmatrix}
   \epsilon & \Delta \\
   \Delta & -\epsilon
   \end{pmatrix},

with bias :math:`\epsilon` and coupling :math:`\Delta` is interacting with a bath of harmonic oscillators defined with,

Bath Hamiltonian:

.. math::

   H_b = \sum_j \left( \frac{P_j^2}{2 M_j}
   + \frac{1}{2} M_j \omega_j^2 R_j^2 \right).

Two-level system bilaterally interacting with a harmonic bath,  defined as System–bath coupling Hamiltonian:

.. math::

   H_{sb} = \sigma_z \sum_j c_j R_j =
   \begin{pmatrix}
   \sum_j c_j R_j & 0 \\
   0 & -\sum_j c_j R_j
   \end{pmatrix}.

Spectral Density
----------------

General definition for :math:`N` bath modes :cite:`Miller1999-debye, Templaar2018`:

.. math::

   J(\omega) = \frac{\pi}{2} \sum_{j=1}^{N}
   \frac{c_j^2}{M_j \omega_j} \, \delta(\omega - \omega_j),

where

  - :math:`M_j` = mass of :math:`j`-th oscillator,  
  - :math:`c_j` = coupling constant between the system and the :math:`j`,-th oscillator,  
  - :math:`\omega_j` = oscillator frequency.

Debye form
^^^^^^^^^^

.. math::

   J_D(\omega) = \frac{E_r}{2} \,
   \frac{\omega \, \omega_c}{\omega^2 + \omega_c^2},

where

  - :math:`E_r` = reorganization energy,  
  - :math:`\omega_c` = characteristic bath frequency.

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

with :math:`\omega_{\text{max}} = 20 \omega_c`, a high-frequency cutoff.

Ohmic form
^^^^^^^^^^

.. math::

   J_O(\omega) = \frac{\pi}{2} \, \alpha \, \omega \, e^{\omega/\omega_c},

where

  - :math:`\alpha` = Kondo parameter that controls the strength of the coupling,  
  - :math:`\omega_c` = characteristic bath frequency.

**Ohmic Bath Discretization**

Frequencies and couplings are discretizeda as:

.. math::

   \omega_j = \omega_c \cdot
   \log\left(1 - \frac{j}{1+N} \right),

.. math::

   c_j = \omega_j \sqrt{
   \frac{M_j \alpha}{1+N}},

