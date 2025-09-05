Tutorials
=========

This section provides step-by-step tutorials for running **SHARP Pack**.  
Each tutorial includes the required ``param.in`` input template, execution commands, and simple plotting instructions.  

.. note::

   - All simulations assume SHARP Pack is compiled (see :ref:`installation`) and the executable is available as ``sharp.x``.  
   - By default, the program reads input parameters from ``param.in`` in the working directory.  
   - Units: All input parameters are in **atomic units**, except temperature in Kelvin. Spin-boson model parameters are given in :math:`\text{cm}^{-1}` when ``enu`` :math:`\neq` ``1``.


**List of tutorials:**

- :ref:`ex1-tully1s`
- :ref:`ex2-tully2p`
- :ref:`ex3-tully3b`
- :ref:`ex4-lchain`
- :ref:`ex5-sb`

.. note::

   All the example tutorials can be found in the ``example/`` directory.

.. note::

   Job submission bash scripts are located in the ``utility/`` directory for the HPC Slurm system, as well as on the local machine. 
