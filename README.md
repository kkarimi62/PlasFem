# PlasFem - Finite Elements Package
Finite element package for modeling amorphous plasticity

## Table of Contents
* [Author Information](#author-information)
* [Compilation](#compilation)
* [Overview](#overview)

## Author Information
* Developer: Kamran Karimi
* Institution: University of Calgary, Department of Physics and Astronomy
* Email: [kamran.karimi1@ucalgary.ca](mailto:lunde@adobe.com?subject=[GitHub]%20Source%20Han%20Sans)

## Compilation
A simpe MakeFile is provided. You just need to run \
\
make \
\
in the package library and the executable file "fem_run" will be created in ./src. To clean old object files (\*.o) and dependencies (\*.d), simply run \
\
make clean \
\
in the library. To delete the executables, run \
\
make clean_all \
\
.

# Overview
The package solves elastic and/or plastic problems by discretizing the domain using two-dimensional linear triangular elements in space. The time
integration is carried out through an explicit scheme besed on the velocity Verlet integrator.
