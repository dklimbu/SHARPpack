.. _installation:

Installation
============

The source code is compiled using the provided ``Makefile`` in the ``source`` directory. The package has been tested under **Linux** with the **Intel Fortran compiler**.

Required Libraries
------------------

To compile and run SHARP Pack, the following external libraries must be installed: **BLAS**, **LAPACK**, and **FFTW3**.

Compilation Steps
-----------------

- Clone the **SHARP Pack 2.0** source code from the  GitHub repository:

  .. code-block:: console

      $ git clone https://github.com/fashakib/SHARP_pack_2.git

- Once the package has been downloaded, navigate to the directory containing the source codes:

  .. code-block:: console

      $ cd SHARP_pack_2/source

- Open the ``Makefile`` and update the library paths to point to your installed libraries.
  For example, edit the lines:

  .. code-block:: make

     LIBPATH = -L/path/to/blas \
             -L/path/to/lapack \
             -L/path/to/fftw3

  .. note::

    Make sure to replace ``/path/to/...`` with the correct installation paths on your system.

- Compile the code by running:

  .. code-block:: console

      $ make

- Upon successful compilation, the ``bin/`` directory will contain the main executable, ``sharp.x``, and the auxiliary executables, ``average.x``, and ``average2db.x``. To verify, run:

  .. code-block:: console

     $ ls ../bin/sharp.x


- Add the ``bin`` directory to your ``PATH`` for convenience:

  .. code-block:: bash

      export PATH=$PATH:/path/to/SHARP_pack_2/bin

Now the **SHARP Pack** is ready to run RPSH simulations.


.. attention::

   When using SHARP Pack, please cite the following papers:

   - D. K. Limbu and F. A. Shakib, `Software Impacts 19, 100604 (2024). <https://doi.org/10.1016/j.simpa.2023.100604>`_
   - D. K. Limbu and F. A. Shakib, `J. Phys. Chem. Lett. 14, 8658–8666 (2023). <https://doi.org/10.1021/acs.jpclett.3c02085>`_
                                                  
