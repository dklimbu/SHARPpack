.. SHARPPack documentation master file, created by
   sphinx-quickstart on Mon Dec  2 08:02:34 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
============

Non-adiabatic dynamics—quantum transitions between electronic and/or vibronic states—is a critical process in electrochemical and photoelectrochemical reactions, including those involved in natural and artificial energy conversion systems. :cite:`Huynh:2007,Hammes-Schiffer:2007`

Simulating these processes accurately in condensed phases with thousands of degrees of freedom (DOFs), especially when nuclear quantum effects (NQEs) like proton tunneling and zero-point energy are significant, is computationally demanding. Mixed quantum-classical dynamics (MQCD) methods address this by treating only light particles quantum-mechanically and treating the environment classically.

One of the most widely used MQCD methods is Fewest-Switches Surface Hopping (FSSH) :cite:`Tully:1990`, where the quantum wavepacket is represented by classical trajectories traveling on a single adiabatic surface interrupted by hops to other surfaces. However, independent classical trajectories used in FSSH are not able to capture decoherence or NQEs accurately.

Approaches like Centroid Molecular Dynamics (CMD) :cite:`Cao:1994` and Ring Polymer Molecular Dynamics (RPMD) :cite:`Craig:2004` offer a way to model NQEs, but have limitations for multi-electron transfer or real-time electronic coherence effects. Ring Polymer Surface Hopping (RPSH) :cite:`Shushkov:2012,Shakib:2017` is a method that combines RPMD and FSSH to overcome these challenges.

SHARP Pack :cite:`Limbu:2024` —short for *Surface Hopping And Ring Polymer* package—is a modular, parallelized, and user-friendly software that implements RPSH, enabling simulations from model systems to more realistic applications.

**Features of SHARP Pack:**

* Supports both FSSH and RPSH simulations.  
* Provides a variety of model Hamiltonians including 1D Tully to spin-boson model with Debye/Ohmic spectral density.
* Handles flexible initial conditions and user-defined control parameters.  
* Implements different schemes for velocity rescaling within the surface hopping algorithm.  
* Propagates dynamics using both velocity-Verlet and Langevin dynamics schemes.  
* Enables serial as well as parallel simulations with trajectory distribution across multiple cores.  
* Offers PIMD sampling with PILE thermostat for initial sampling.
* Option for Diagonal Born-Oppenheimer correction (DBOC) for accurate quantum dynamics.
* Presents option for decoherence correction to improve population and coherence evolution.
* Produces a wide range of outputs: energies, NACs, populations, branching probablities, etc.  
* Provides ready-to-use job scripts for local runs, HPC cluster submissions, and branching probability analysis.  
* Supplies comprehensive tutorials and example cases with scripts for visualization.  

