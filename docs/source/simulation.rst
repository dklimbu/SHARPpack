.. SHARPPack documentation master file, created by
   sphinx-quickstart on Mon Dec  2 08:02:34 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Simulations
===========

To run a simulation with **SHARP Pack**, only a single input file named
``param.in`` is required. This file contains all simulation details specified through predefined keywords.

.. note::
   The file must be located in the same directory where the executable
   ``sharp.x`` is launched, unless a path is explicitly provided.

Keywords
--------

Table :ref:`tab-keyword` lists all available input parameter keywords and their usage.

.. _tab-keyword:

.. list-table:: Input Keywords
   :widths: 20 15 65
   :header-rows: 1
   :align: center

   * - **Keywords**
     - **Arguments**
     - **Descriptions**
   * - **acckval**
     - `yes`
     - If TRUE, writes accumulated result with k-value; useful for branching probability calculation(s) of **Tully model**
   * -
     - :math:`\color{blue}{\text{no}}`
     - No accumulated result(s)
   * - **approximation**
     - :math:`\color{blue}{\text{CA}}`
     - Centroid approximation RPSH (RPSH-CA)
   * - 
     - BA
     - Bead approximation RPSH (RPSH-BA)
   * - **dboc**
     - :math:`\color{blue}{\text{BF}}`
     - DBOC corrected force(s) at bead level
   * - 
     - CF
     - DBOC corrected force(s) at centroid level
   * - **decoherence**
     - :math:`\color{blue}{\text{none}}`
     - No decoherence correction
   * - 
     - damp
     - Expontial damping decoherence 
   * - **dynamics**
     - :math:`\color{blue}{\text{rpmd}}`
     - RPMD dynamics with surface hopping
   * - 
     - :math:`\text{trpmd}`
     - T-RPMD dynamics with surface hopping
   * - **iprint**
     - :math:`n`
     - Print data every n timesteps
   * - **iseed**
     - :math:`n`
     - Random number generator seed
   * - **nbead**
     - :math:`n`
     - Number of bead(s)
   * - **ncore**
     - :math:`n`
     - Number of core(s): (1-Serial), (n-Parallel)
   * - **nequil**
     - :math:`n`
     - equailibrate steps **if pimd sampling**
   * - **nmode**
     - :math:`\color{blue}{\text{fft}}`
     - FFT tranform for normal mode
   * -
     - matrix 
     - use normal mode tranformation matrix for normal mode
   * - **nparticle**
     - :math:`n`
     - Number of simulating particle(s)
   * - **nsample**
     - :math:`n`
     - sample run steps **if pimd sampling**
   * - **nstate**
     - :math:`\color{blue}{n}`
     - Number of electronic states; Default based on Model
   * - **nstep**
     - :math:`n`
     - Number of simulation steps
   * - **ntraj**
     - :math:`n`
     - Number of trajectories
   * - **model**
     - :math:`s`
     - Model name string from Table :ref:`tab-model`
   * - **pimd**
     - :math:`\text{pile} \quad \tau`
     - PIMD sampling of postion and velocity with PILE thermostat
   * - **pinit**
     - :math:`f`
     - Initial momentum (a.u.)
   * - **rsamp**
     - fraction
     - Sample centroid bead from Gaussian distribution and other  bead(s) position based on fraction of de-Broglie length
   * -
     - :math:`\color{blue}{\text{gaussian}}`
     - Sample initial bead(s) position from Gaussian distribution
   * -
     - wigner
     - Sample initial bead(s) position from Wigner distribution
   * - **rmap**
     - :math:`\color{blue}{\text{no}}`
     - No mapping ring polymer centroid to dividing surface after sampling **rsamp**
   * - 
     - r0
     - Mapping ring polymer centroid to dividing surface after sampling **rsamp**
   * - **rmode**
     - :math:`\color{blue}{\text{direct}}`
     - Direct ring polymer sampling of **rsamp** 
   * - 
     - norm
     - Normal mode ring polymer sampling of **rsamp** for **gaussian**
   * - **rundtail**
     - `yes`
     - Print run-time details of the simulation
   * -
     - :math:`\color{blue}{\text{no}}`
     - No details of simulation
   * - **temperature**
     - :math:`f`
     - Temperature (K)
   * - **tstep**
     - :math:`f`
     - Simulation time-step (a.u.)
   * - **vsamp**
     - fixed
     - Deterministic initial momentum
   * -
     - :math:`\color{blue}{\text{gaussian}}`
     - Sample initial momenta from Gaussian distribution
   * -
     - wigner
     - Sample initial momenta from Wigner distribution
   * - **vrescale**
     - :math:`\color{blue}{\text{CL}}`
     - Centroid-approximation level of velocity rescaling of ring-polymer for energy conservation
   * -
     - BL
     - Bead-approximation level of velocity rescaling of ring-polymer for energy conservation
   * - **vreverse**
     - :math:`\color{blue}{\text{never}}`
     - Never reversal (NR) velocity for frustrated hop(s)
   * -
     - always
     - Always reversal (VR) of velocity for frustrated hop(s)
   * -
     - delV1
     - Velocity reversal based on Truhlar’s :math:`\Delta V` scheme
   * -
     - delV2
     - Velocity reversal based on Subotnik’s :math:`\Delta V^2` scheme
   * - **finish**
     - `***` 
     - End of the input file

..   * - **usrkval**
     - yes
     - If TRUE, read initial momentum (k) value and overwrite `param.in` k-value
   * -
     - :math:`\color{blue}{\text{no}}`
     - Read and use initial k-value from `param.in`

*Legend:* 
   :math:`n` – integer, 
   :math:`f` – real number, 
   :math:`s` – string,
   :math:`\color{blue}{\text{blue}}` – denotes default value.

----

.. list-table:: Spin-boson model specific keywords
   :widths: 20 30 65
   :header-rows: 1
   :align: center

   * - **Keywords**
     - **Arguments**
     - **Descriptions**
   * - **model**
     - :math:`\text{spinboson} \quad \epsilon \quad \Delta \quad \text{enu}`
     - Spin-boson model with parameters with scale factor ``enu``
   * - **spectra**
     - debye :math:`\quad E_r \quad \omega_c \quad T`
     - Debye spectral density parameters and temperature
   * - 
     - ohmic :math:`\quad E_r \quad \omega_c \quad T`
     - Ohmic spectral density parameters and temperature

.. important::

   - All input parameters are specified in **atomic units (a.u.)**, except for **temperature**, which must be given in **Kelvin**.  

   - For **spin-boson parameters**, values are expressed in energy unit of :math:`\text{cm}^{-1}` including temperature, and scaled by the factor ``enu``.  

   - If ``enu = 1``, the spin-boson parameters are interpreted directly as in **atomic units**.


Input File
----------

A sample input file (param.in) for Spin-Boson model.

   .. literalinclude:: param.in
      :language: bash
      :caption: ``SAMPLE INPUT FILE (param.in)``

Output Files
------------

After a successful simulation run, SHARP Pack generates the following output files:

.. _tab-outfile:

.. list-table:: Output file(s)
   :widths: 30 100
   :header-rows: 1
   :align: center

   * - **File(s)**
     - **Descriptions**
   * - ``param.out``
     - Lists all the model parameters used in the simulation.
   * - ``dcoupling.out``
     - Detailed trajectory data including positions, momenta, energies, NACs, etc.
   * - ``hoppinghist.out``  
     - Records surface-hopping events between different potential energy surfaces.
   * - ``pop_adiabat1.out``  
     - Adiabatic population computed using Method I.
   * - ``pop_adiabat2.out``  
     - Adiabatic population computed using Method II.
   * - ``pop_diabat1.out``  
     - Diabatic population computed using Method I.
   * - ``pop_diabat2.out``  
     - Diabatic population computed using Method II.
   * - ``pop_diabat3.out``  
     - Diabatic population computed using Method III.
   * - ``pop_branch.out``  
     -  Branching probablity of **Tully model** if ``accKval``  ``yes`` .
   * - ``sampling_pos.out``  
     - Initial sampling distribution of positions.  
   * - ``sampling_vel.out``  
     - Initial sampling distribution of velocities.  
   * - ``sampling_dis.out``  
     - Initial sampling Boltzmann distribution function of speeds.  
   * - ``pimd_sample.out``  
     - PIMD sampling output (generated only if PIMD sampling is enabled).  
   * - ``energy_surface.out``  
     - Diabatic and adiabatic energy surfaces, including NAC vectors.  
   * - ``energy_dboc.out``  
     - Adiabatic and DBOC energy surfaces (if dboc is requested)
   * - ``psi.out``  
     - Wavefunctions, adiabatic coefficients, and PES.  
   * - ``bathfrequency.out``  
     - Spectral density information (specific to the spin-boson model).  

.. note::

   The output files ``dcoupling.out``, ``hoppinghist.out``, and ``psi.out`` are produced **only if** ``rundtail`` is ``yes`` in the input file.
   

Population Calculation
----------------------

SHARP Pack runs simulations in an **adiabatic formalism**; hence, the adiabatic state populations are the direct results. But it also provides multiple approaches for computing adiabatic and diabatic populations.  

**Adiabatic Population**
^^^^^^^^^^^^^^^^^^^^^^^^

There are two methods:  

 **Method-I: Direct counting trajectory**

 The percentage of trajectories propagating on each adiabatic PES represents the population of that state.

   .. math::

      P_{\alpha} = \langle \delta_{\alpha,\lambda} \rangle 

   where :math:`\lambda` is the adiabatic active surface.  

   output file name: ``pop_adiabat1.out``  

 **Method-II: Density matrix method**

 Populations are computed from the complex-valued electronic coefficients.

   .. math::

      P_{\alpha} = \langle \rho_{\alpha\alpha} \rangle

   where :math:`\rho_{ij} = c_i\,c_j^*` is the electronic density matrix.  

   output file name: ``pop_adiabat2.out``  


**Diabatic Population**
^^^^^^^^^^^^^^^^^^^^^^^

There are three methods, :cite:`Subotnik2013:populaiton` using :math:`U_{ij}` which are elements of the adiabatic-to-diabatic transformation matrix, and :math:`\lambda` denotes the active surface.

 **Method-I: Projection Method**

 The diabatic population is determined by projecting the active adiabatic state onto the diabatic basis.

   .. math::

      P_\alpha = \left\langle \sum_i |U_{\alpha i}\, \delta_{i,\lambda}|^2 \right\rangle,

   output file name: ``pop_diabat1.out``  

 **Method-II: Electronic Wavefunction Method**
 
 Diabatic populations are calculated by transforming the electronic coefficients :math:`c_i` in the adiabatic basis to the diabatic representation.

   .. math::

      P_\alpha = \left\langle \sum_i | U_{\alpha i}\,c_i |^2 \right\rangle,

   output file name: ``pop_diabat2.out``  


 **Method-III: Mixed quantum-classical density matrix**
 
 The diabatic populations are extracted from the electronic density matrix in the diabatic representation.

   .. math::

      P_\alpha = \left\langle \sum_i |U_{\alpha i}|^2 \delta_{i,\lambda} + \sum_{i<j} 2\, \text{Re}\left[ U_{\alpha i} , \rho_{ij} , U_{\alpha j}^* \right] \right\rangle,

   output file name: ``pop_diabat3.out``  

Running Simulation
------------------

To run **SHARP Pack**, a single input file ``param.in`` is required, which defines the system model parameters (see :ref:`tab-keyword` for details). Simulations are executed using the provided bash scripts located in the ``utility/`` directory.

Available job submission scripts:

- ``job-script-local.sh``  
  Run simulations on a **local machine**.

- ``job-script-hpc.sh``  
  Submit jobs on an **HPC cluster** using the **Slurm** workload manager.

- ``job-script-branching.sh``  
  Used **only** for computing branching probabilities on **Tully model**  over multiple *k*-values on an HPC system.

- ``job-script-avg.sh`` 
  Submit job that computes average of `output` of the completed parallel jobs simulation.

Each script automatically generates the appropriate job submission file (serial or parallel, depending on the ``ncore`` value in ``param.in``) and submits the job.

**Usage**

Run any of the scripts as follows:

.. code-block:: bash

   # Run on a local machine
   $ sh job-script-local.sh

   # Run on an HPC cluster (Slurm)
   $ sh job-script-hpc.sh

   # Compute branching probability on HPC
   $ sh job-script-branching.sh

----

**Example script**

An example of ``job_script_hpc.sh`` is shown below which is designed for running simulations on an **HPC cluster** using the Slurm workload manager.
It automatically detects the number of cores (``ncore``) requested in the ``param.in`` file and prepares the submission accordingly:

- If ``ncore  1``; A **serial job** is created and submitted.

- If ``ncore  n`` (``n > 1``); The script creates n **run* directories** and distributes the number of trajectories across them for **parallel simulation**.
  After all simulations are complete, the results are averaged to produce the final output.

----

.. literalinclude:: job_script.sh
    :language: bash
    :caption: ``Bash script for job submission``

.. note::

   To compute the average results from parallel runs, must to specify ``fname``, and ``df`` inside the **write_average_input()** fucntion in ``job_script_hcp.sh`` or ``job_script_local.sh`` as:  

   .. code-block:: bash

       # ==== FUNCTION TO CREATE AVERAGE INPUT ====
       write_average_input() {
       #fname :: File name of specific population data for averaging
       #df    :: No. of data column(s)(except first column)

       fname="pop_diabat3.out"
       df=2
       .
       .

   The final average output file is ``fname_ave.out``. e.g. ``pop_diabat3.out`` ==> ``pop_diabat3_ave.out``

