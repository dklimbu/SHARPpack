.. _changelog:

Change Log
==========

Notable changes to **SHAPRP Pack** will be documented here.


Version 2.0 (2025-xx-xx)
------------------------

New Features
^^^^^^^^^^^^

* Extended method to bead-approximation ring polymer surface hopping (RPSH-BA).
* Diagonal Born-Oppenheimer Correction (DBOC) implementation. 
* Option for decoherence correction.
* Random seed handling updated for reproducibility.
* Added PILE thermostat for PIMD sampling simulations.
* New PIMD and runtime outputs: `bead`- and `centroid`-level system temperatures.  
* New option for constrained centroid to dividing surface (``R=0``) for initial sampling, shifting the whole ring polymer accordingly.  
* Optput histogram distribution of initial sampling positions and momenta.  
* Extented model system to *spin-boson* and *dboc1* models.
* Added *Debye/Ohmic* spectral density options for *spin-boson model*.  
 

Bug Fixes
^^^^^^^^^

* Initial bead position sampling from de Broglie length for SB model
* Initial force calculation at t=0 carried out before running dynamics.  

Documentation
^^^^^^^^^^^^^

- Online documentation lauch (docs in Sphinx via Readthedocs)

Version 1.0 (2023-11-28)
------------------------

New Features
^^^^^^^^^^^^

- SHARP Pack software introduce based on Ring Polymer Surface Hopping method.

