.. SHARPPack documentation master file, created by
   sphinx-quickstart on Mon Dec  2 08:02:34 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SHARP Pack documentation
========================

.. figure:: _static/rpsh.png
   :scale: 20 %
   :align: center
   :alt: SHARPPack

Welcome to **SHARP Pack** documentation, 

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. note::
   This website is in under active development...

Introduction
============

Non-adiabatic dynamics—quantum transitions between electronic and/or vibronic states—is a critical process in electrochemical and photoelectrochemical reactions, including those involved in natural and artificial energy conversion systems.

Simulating these processes accurately in condensed phases with thousands of degrees of freedom (DOFs), especially when nuclear quantum effects (NQEs) like proton tunneling and zero-point energy are significant, is computationally demanding.

**Mixed quantum-classical dynamics (MQCD)** methods address this by treating only light particles quantum-mechanically and treating the environment classically.

One of the most widely used MQCD methods is **Fewest-Switches Surface Hopping (FSSH)**, where the quantum wavepacket is represented by classical trajectories traveling on a single adiabatic surface interrupted by hops to other surfaces.

However, independent classical trajectories used in FSSH are not able to capture decoherence or NQEs accurately.

Approaches like **Centroid Molecular Dynamics (CMD)** and **Ring Polymer Molecular Dynamics (RPMD)** offer a way to model NQEs, but have limitations for multi-electron transfer or real-time electronic coherence effects.

**Ring Polymer Surface Hopping (RPSH)** is a method that combines RPMD and FSSH to overcome these challenges.

SHARP Pack—short for *Surface Hopping And Ring Polymer* package—is a modular, parallelized, and user-friendly software that implements RPSH, enabling simulations from model systems to more realistic applications.

Features of SHARP Pack:

- Coupled electronic-nuclear motion
- Surface hopping algorithm
- Feynman's path integral formalism
- Extensible for new methodologies and models
- Efficient yet accurate nonadiabatic dynamics

Models
======

.. _tab-model:

.. list-table:: Models implemented in SHARP Pack
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
   * - **Morse Models**
     -
   * - `morse1`
     - Morse Model I – photo-dissociation
   * - `morse2`
     - Morse Model II – photo-dissociation
   * - `morse3`
     - Morse Model III – photo-dissociation
   * - **Detailed Balance Models**
     -
   * - `db2lchain`
     - 2-state model coupled with linear chain
   * - `db3lchain`
     - 3-state super-exchange model with linear chain coupling

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
^^^^^^^^^^^^^^^^^^^^^^

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


Spin-Boson Model
^^^^^^^^^^^^^^^^

Total Hamiltonian:

.. math::

   H = H_s + H_b + H_{sb}

Where:

- :math:`H_s`: Spin system (2-level)
- :math:`H_b`: Harmonic bath
- :math:`H_{sb}`: System-bath coupling

Spectral density:

.. math::

   J(\omega) = \frac{E_r}{2} \cdot \frac{\omega \omega_c}{\omega^2 + \omega_c^2}

Discretized mode sampling:

.. math::

   \omega_j = \tan\left(\frac{j}{N} \tan^{-1}\left(\frac{\omega_{\text{max}}}{\omega_c}\right)\right) \cdot \omega_c

.. math::

   c_j = \omega_j \sqrt{\frac{E_r}{\pi N} \tan^{-1}\left(\frac{\omega_{\text{max}}}{\omega_c}\right)}

Sample Initialization
---------------------

Tully Model
^^^^^^^^^^^

Gaussian wavepacket:

.. math::

   \langle R|\chi(R)\rangle = \left(\frac{2\alpha}{\pi}\right)^{1/4} \exp\left[-\alpha(R - R_0)^2 - \frac{iR(P - P_0)}{\hbar} \right]

Wigner distribution:

.. math::

   \rho^W_R = \frac{1}{\pi \hbar} \exp\left[-2\alpha(R - R_0)^2 - \frac{(P - P_0)^2}{2\alpha \hbar^2} \right]

Constants:

.. list-table:: Constants for Tully Models
   :widths: 30 70
   :header-rows: 1

   * - **Parameter**
     - **Value**
   * - :math:`\alpha`
     - 0.25
   * - :math:`m`
     - 2000
   * - :math:`\sigma_R`
     - :math:`\sqrt{1 / 2\alpha}`
   * - :math:`\sigma_P`
     - :math:`\sqrt{\hbar^2 \alpha / 2}`

Morse Model
^^^^^^^^^^^

.. list-table:: Morse Model Initialization
   :widths: 30 70
   :header-rows: 1

   * - **Parameter**
     - **Value**
   * - :math:`m`
     - 20000
   * - :math:`\omega_0`
     - 0.005
   * - :math:`\sigma_R`
     - :math:`\sqrt{\frac{\hbar}{2 m \omega_0}}`
   * - :math:`\sigma_P`
     - :math:`\sqrt{\frac{m \hbar \omega_0}{2}}`

Spin-Boson Model
^^^^^^^^^^^^^^^^

Boltzmann sampling:

.. math::

   P(p, q) \propto \prod \exp\left[-\beta \left(\frac{p^2}{2m} + \frac{1}{2} m \omega^2 q^2 \right)\right]

Wigner sampling:

.. math::

   \sigma^w_R = \sqrt{\frac{\hbar}{2m\omega \tanh(\beta \hbar \omega / 2)}}

.. math::

   \sigma^w_P = \sqrt{\frac{m \hbar \omega}{2 \tanh(\beta \hbar \omega / 2)}}

Frustrated Hops in Tully Model
------------------------------

Velocity rescaling after hopping:

.. math::

   \dot{R}'_i = \dot{R}_i - \gamma \frac{d_{kk'}^i}{m_i}

Velocity reversal strategies:

.. list-table:: Velocity Reversal Options
   :widths: 30 70
   :header-rows: 1

   * - **Option**
     - **Description**
   * - `never`
     - Never reverse momentum
   * - `always`
     - Always reverse on failed hop
   * - `truhlar`
     - Use Truhlar's criterion
   * - `subotnik`
     - Use Subotnik's criterion

Compilation
-----------

Requirements:

- Linux OS
- Intel Fortran compiler
- BLAS, LAPACK, FFTW3

Build:

.. code-block:: bash

   cd source/
   make

Output: `bin/sharp.x`


.. list-table:: Input Keywords
   :widths: 20 15 65
   :header-rows: 1
   :align: center

   * - **Keywords**
     - **Arguments**
     - **Descriptions**
   * - **model**
     - :math:`s`
     - Model name string from Table :ref:`tab-model`
   * - **nbeads**
     - :math:`n`
     - Number of bead(s)
   * - **nparticle**
     - :math:`n`
     - Number of simulating particle(s)
   * - **nsteps**
     - :math:`n`
     - Number of simulation steps
   * - **ntraj**
     - :math:`n`
     - Number of trajectories
   * - **nstates**
     - :math:`\color{blue}{n}`
     - Number of electronic states; Default based on Model
   * - **ncpu**
     - :math:`n`
     - Number of CPU(s) (1-Serial), (n-Parallel)
   * - **timestep**
     - :math:`f`
     - Simulation time-step (a.u.)
   * - **temp**
     - :math:`f`
     - Temperature (K)
   * - **kinitial**
     - :math:`f`
     - Initial momentum (a.u.)
   * - **beadpos**
     - `same`
     - Assign same position for bead(s) of each particle
   * -
     - :math:`\color{blue}{\text{gaussian}}`
     - Initial (beads) position from Gaussian distribution
   * -
     - `wigner`
     - Initial (beads) position from Wigner distribution
   * - **vinitial**
     - `fixed`
     - Deterministic initial momentum
   * -
     - :math:`\color{blue}{\text{gaussian}}`
     - Initial momentum from Gaussian distribution
   * -
     - `wigner`
     - Initial momentum from Wigner distribution
   * - **vreverse**
     - :math:`\color{blue}{\text{never}}`
     - Never reversal (NR) velocity for frustrated hop(s)
   * -
     - `always`
     - Always reversal (VR) of velocity for frustrated hop(s)
   * -
     - `truhlar`
     - Velocity reversal based on Truhlar’s :math:`\Delta V` scheme
   * -
     - `subotnik`
     - Velocity reversal based on Subotnik’s :math:`\Delta V^2` scheme
   * - **vrescale**
     - :math:`\color{blue}{\text{ca}}`
     - Centroid-approximation level of velocity rescaling of ring-polymer for energy conservation
   * -
     - `ba`
     - Bead-approximation level of velocity rescaling of ring-polymer for energy conservation
   * - **usrkval**
     - `yes`
     - If TRUE, read initial momentum (k) value and overwrite `param.in` k-value
   * -
     - :math:`\color{blue}{\text{no}}`
     - Read and use initial k-value from `param.in`
   * - **acckval**
     - `yes`
     - If TRUE, writes accumulated result with k-value; useful for branching probability calculation(s)
   * -
     - :math:`\color{blue}{\text{no}}`
     - No accumulated result(s)
   * - **rundtail**
     - `yes`
     - Print run-time details of the simulation
   * -
     - :math:`\color{blue}{\text{no}}`
     - No details of simulation
   * - **iprint**
     - :math:`n`
     - Print data every n timesteps
   * - **finish**
     - *(empty)*
     - End of the input file

*Legend:* :math:`n` – integer, :math:`f` – real number, :math:`s` – string,
:math:`\color{blue}{\text{blue}}` – denotes default value.



Input File: param.in
--------------------

Sample configuration:

.. code-block:: text

   model       tully1
   method      rpsh 4
   nbeads      4
   ntraj       100
   temp        300
   timestep    0.1
   vreverse    subotnik
   finish

Output Files
------------

.. list-table:: Output Files
   :widths: 30 70
   :header-rows: 1

   * - **File**
     - **Description**
   * - `param.out`
     - All model parameters
   * - `dcoupling.out`
     - Simulation log
   * - `hopping.out`
     - Surface hopping events
   * - `adiabat[1-2].out`
     - Adiabatic populations
   * - `diabat[1-3].out`
     - Diabatic populations
   * - `pkval.out`
     - Branching probability (if enabled)

Running Simulations
-------------------

.. code-block:: bash

   ./sharp.x

or

.. code-block:: bash

   sh submit.sh

Examples
--------

**Example 1: Tully I (serial)**

.. code-block:: bash

   cd example/tully1/
   sh submit.sh

**Example 2: Tully II (parallel)**

Edit `param.in`:

.. code-block:: text

   model     tully2
   method    rpsh 4
   ncpu      10

Run:

.. code-block:: bash

   sh submit.sh

**Example 3: Tully III with branching probability**

.. code-block:: text

   model     tully3
   usrkval   yes
   acckval   yes

Run:

.. code-block:: bash

   sh job.sh

**Example 4: db2lchain with detailed balance**

.. code-block:: text

   model     db2lchain
   method    fssh
   ncpu      20
   temp      1000

Plot results:

.. code-block:: bash

   gnuplot plot4.in




