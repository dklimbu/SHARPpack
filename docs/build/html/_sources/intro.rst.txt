.. SHARPPack documentation master file, created by
   sphinx-quickstart on Mon Dec  2 08:02:34 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
============

.. figure:: _static/rpsh.png
   :scale: 20 %
   :align: center
   :alt: SHARPPack

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
