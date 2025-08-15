.. SHARPPack documentation master file, created by
   sphinx-quickstart on Mon Dec  2 08:02:34 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Keywords
========

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
   * - **nbead**
     - :math:`n`
     - Number of bead(s)
   * - **nparticle**
     - :math:`n`
     - Number of simulating particle(s)
   * - **nstep**
     - :math:`n`
     - Number of simulation steps
   * - **ntraj**
     - :math:`n`
     - Number of trajectories
   * - **nstate**
     - :math:`\color{blue}{n}`
     - Number of electronic states; Default based on Model
   * - **ncpu**
     - :math:`n`
     - Number of CPU(s) (1-Serial), (n-Parallel)
   * - **tstep**
     - :math:`f`
     - Simulation time-step (a.u.)
   * - **temperature**
     - :math:`f`
     - Temperature (K)
   * - **kinitial**
     - :math:`f`
     - Initial momentum (a.u.)
   * - **beadpos**
     - same
     - Assign same position for bead(s) of each particle
   * -
     - :math:`\color{blue}{\text{gaussian}}`
     - Initial (beads) position from Gaussian distribution
   * -
     - wigner
     - Initial (beads) position from Wigner distribution
   * - **vinitial**
     - fixed
     - Deterministic initial momentum
   * -
     - :math:`\color{blue}{\text{gaussian}}`
     - Initial momentum from Gaussian distribution
   * -
     - wigner
     - Initial momentum from Wigner distribution
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
   * - **vrescale**
     - :math:`\color{blue}{\text{CA}}`
     - Centroid-approximation level of velocity rescaling of ring-polymer for energy conservation
   * -
     - BA
     - Bead-approximation level of velocity rescaling of ring-polymer for energy conservation
   * - **usrkval**
     - yes
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
   * - **iseed**
     - :math:`n`
     - Random number generator seed
   * - **dboc**
     - :math:`\color{blue}{\text{bead}}`
     - DBOC corrected force(s) at bead level
   * - 
     - centroid
     - DBOC corrected force(s) at centroid level
   * - **decoherence**
     - damp
     - Expontial damping decoherence 
   * - 
     - :math:`\color{blue}{\text{none}}`
     - No decoherence correction
   * - **approximation**
     - :math:`\color{blue}{\text{CA}}`
     - Centroid approximation RPSH (RPSH-CA)
   * - 
     - BA
     - Bead approximation RPSH (RPSH-BA)
   * - **nmode**
     - :math:`\color{blue}{\text{fft}}`
     - FFT tranform for normal mode
   * -
     - matrix 
     - use normal mode tranformation matrix for normal mode
   * - **finish**
     - `***` 
     - End of the input file

*Legend:* 
   :math:`n` – integer, 
   :math:`f` – real number, 
   :math:`s` – string,
   :math:`\color{blue}{\text{blue}}` – denotes default value.



.. list-table:: Spin-boson model specific keywords
   :widths: 20 15 65
   :header-rows: 1
   :align: center

   * - **Keywords**
     - **Arguments**
     - **Descriptions**
   * - **model**
     - :math:`\text{spinboson} \quad \epsilon \quad \Delta \quad \text{enu}`
     - Spin-boson model with parameters
   * - **spectra**
     - debye :math:`\quad E_r \quad \omega_c \quad T`
     - Debye spectral density
   * - 
     - ohmic :math:`\quad E_r \quad \omega_c \quad T`
     - Ohmic spectral density
