.. _ex2-tully2p:

Tully Model-II (RPSH)
=====================

This tutorial demonstrates how to run RPSH simulation in **parallel mode**.
The number of parallel runs is controlled by the value of ``ncore  n`` specified in the ``param.in`` file.

Follow the steps below to perform the simulation:

Prepare ``param.in``
--------------------

Set ``ncore  n`` (where ``4`` is the number of desired parallel runs).

  .. code-block:: bash

    #tully model-II input parameters
    model        tully2
    nParticle    1
    nbeads       4
    approx       CA
    nmode        matrix
    iseed        12345
    ncore        4
    tstep        1.0
    nsteps       3000
    ntraj        1000
    pinit        15.0
    rsamp        gaussian
    vsamp        fixed
    rundtail     Yes
    iprint       10
    finish


Run Simulation
--------------

Run the simulation using the provided job script in ``utility/`` directory (e.g., ``job-script-hpc.sh`` for an HPC cluster, or ``job-script-local.sh`` for local execution).  

Based on ``ncore  4``, the job script will:

- Create **4** ``run*`` directories.
- Distribute the trajectories evenly across the runs.  
- Execute the simulations in parallel.
- After all parallel simulations are complete, the script automatically averages the results from the **4** ``run*`` and produces the final population output.


  .. note::
    To get average result from parallel runs, must to specify ``fname``, and ``df`` inside the **write_average_input()** fucntion in ``job_script_hcp.sh`` or ``job_script_local.sh`` as:  

    .. code-block:: bash

       # ==== FUNCTION TO CREATE AVERAGE INPUT ====
       write_average_input() {
       #fname :: File name of specific population data for averaging
       #df    :: No. of data column(s)(except first column)

       fname="pop_diabat3.out"
       df=2
       .
       .

  .. note::
    To compute a standalone average of a specific output file after completing the parallel simulation, one can use ``job_script_avg.sh`` script provided in ``utility/`` directory.  

  .. warning::   
    Pay attention to modifying the variables ``fname`` (output filename) and ``df`` (number of data columns, excluding the first column) inside the ``write_average_input()`` function of ``job_script_avg.sh``.

    .. code-block:: bash

       .
       .
       # ==== FUNCTION TO CREATE AVERAGE INPUT ====
       write_average_input() {
       #fname :: File name of specific population data for averaging
       #df    :: No. of data column(s)(except first column)

       fname="pop_diabat2.out"
       df=4
       .
       .

       

Plot Results
------------

Finally, use any visualization tool (e.g., **gnuplot** or **matplotlib**) to plot the population dynamics from the output files and analyze the results.  

  .. code-block:: gnuplot

    #!/usr/bin/gnuplot

    set xlabel "Time (a.u.)"
    set ylabel "Population"
    set title "Tully2 Diabatic Populations"
    plot "pop_diabat3_ave.out" u 1:2 w l lw 2 title "Diabatic State 1", \
         "pop_diabat3_ave.out" u 1:3 w l lw 2 title "Diabatic State 2"

