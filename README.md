# SHARP pack software package
**SHARP pack**, short for _**Surface Hopping And Ring Polymer package**_, is a highly-parallelized software with a modular structure that is developed to (1) help the development of approximate quantum dynamics methods and (2) simulate non-adiabatic dynamics in condensed phases subject to nuclear quantum effects.

## Installation
The SHARP pack software package can be downloaded from the Github repository. Create a folder and git clone this repository.
```
$ git clone https://github.com/dklimbu/SHARPpack.git
```

Once the package has been downloaded, navigate to the directory containing the package and run the following command:
```
$ make
```

## Prerequisite
The SHARP pack is tested under a Linux environment with an Intel Fortran compiler. Path for BLAS, LAPACK, and FFTW3 libraries need to be provided.

## Usage
The SHARP pack software package includes several modules that can be used to perform different types of simulations and analyses. See the SHARP Pack manual for more detail (coming soon).

## Examples
The SHARP pack software package includes several examples that demonstrate how to use the different method(s)/model(s). These examples can be found in the example/ directory of the package.

## Running Simulation
After compiling, **sharp.x** can be called upon to run the simulation in a serial mode
(**ncpu 1**). An input file, _**‘param.in’**_, is required to run any SHARP pack simulation.
Alternatively, the SHARP pack can be run using the job submission bash script (see
_**submist.sh**_ script in **bin/** directory). Based on ncpu in _**‘param.in’**_, the code will run as _serial_
on a single node or _parallel_ (trajectories are parallelized in this case) on nodes on the cluster.

## Contact
For any queries and feedback please contact Dr. Limbu (dil.limbu@umkc.edu) or Dr. Shakib (shakib@njit.edu).