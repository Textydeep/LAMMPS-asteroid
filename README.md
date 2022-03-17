# LAMMPS-asteroid
A LAMMPS module to run simulations of granular flow over rotating and gravitating tri-axial ellipsoids without tumble.

1. First install the current version of LAMMPS.
2. Install GSL locally on your allocation or have it centrally installed on your super-computer.
3. If locally installed make sure the path name inside fix_ElGRAV.cpp is modified appropriately. You may also need to modify path names in the dependencies of the fix_ELGRAV.cpp and their dependencies.
4. If installed centrally on the cluster, things are much more easy and you may just need to include the gsl headers in code.
5. Makefile may (central installation, will probably need to add flags) or may not (local installation) need modifications. Otherwise, any module with the name of the form 'fix_modulename.cpp' and 'fix__modulename.h' will be compiled automatically by the Makefile.

The module does the following: Computes gravity field of a tri-axial ellipsoid using elliptic integrals from GSL. For that it also solves a cubic. The math can be found here: https://www.deepayanbanik.info/_files/ugd/7a3bd6_02b6d98e8e75443b9762b334450946c6.pdf. 

Additionally, it also solves the particle dynamics in a rotating frame of reference. The rotation pole is fixed and oriented along the z-axis of the simulation box in which you have your asteroid. So if you want to simulate a rotating asteroid with grains flowing on the surface, then fix the core (name these particles as a group, say group 1) using the command 'freeze' in LAMMPS and let the flowing particles (say, group 2) be acted on by the 'ElGRAV' command.

The 'ElGRAV' command is similar to the 'gravity' command but has some minor modifications. A typical syntax would be:

fix    grav 20 ElGRAV 80000 vector 250 250 225 0

where, 
grav ---> name of fix
20  ----> group id of particles which will be subjected to gravity and rotation
80000 --> density of the core
250 ----> semi-x axis of ellipsoid
250 ----> semi-y axis of ellipsoid
225 ----> semi-z axis of ellipsoid
0 ------> spin rate of ellipsoid
