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
The package contains 3 independent types of programs: (i) A program that derives a density matrix from a real earthquake catalogue (Chad's format)
   for an inhomogeneous background-earthquake distribution, (ii) the actual
   ETAS catalogue generator that can use a previously generated density matrix
   OR generate random (homogeneous) background earthquakes, and (iii) a tool to
   clean up the generated catalogues regarding detector blindness after very
   strong earthquakes.
