# SGporecondensation
Small repository containing the relevant simulation and analysis code for stress-granule condensation at pore damaged vesicles (part of "Stress granules plug and stabilize damaged endolysosomal membranes" - DOI: 10.1038/s41586-023-06726-w). For the study of plugging effects.

## Simulation note
We recommend using the LAMMPS version from November 3rd 2022 or a later one. The ylz potential (which is not present in older versions) is required for these simulations. To run a simulation just execute the LAMMPS executable (lmp_mpi / lmp_serial) within the simulation folder as: /$YOURPATHTOLAMMPS/lmp_mpi - in in.local

## Analysis note
For the analysis performed the surface reconstruction library developped by Miguel Amaral (Saric lab) is required. This is freely available at https://github.com/Saric-Group/monolayer_shared
