Cantera-CounterflowFlame
==================

## Description

This package is part of the flameletFoam package and provides a tool to generate the tables needed for the simulations.
The program uses the Cantera chemistry libraries to solve the chemical kinetics for a one-dimensional counterflow diffusion flame.
The output of the solver is a flamelet table that can be read by the OpenFOAM utility canteraToFoam and then used by the CFD solver flameletFoam

More information is available on the Extend-bazaar page:
https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/flameletFoam

## Installation

This version works with Cantera 2.0.0.

* Prepare a directory on your system, e.g.:

  `mkdir ~/cantera/Cantera-CounterflowFlame/`

* Download flameletFoam using git:

  `git clone https://github.com/flameletFoam/Cantera-CounterflowFlame/ ~/cantera/Cantera-CounterflowFlame/`


* Set the correct path to your Cantera installation in ./src/Makefile, e.g replace the line 

  `($YOUR_CANTERA_INSTALLATION_FOLDER$)/build/platform/Cantera.mak`
 
  by

  `~/cantera/cantera-2.0.0/build/platform/Cantera.mak`

* Build the source by executing `make` in the src folder.

If everything worked, there should be an executable `flamelet` in the parent folder.

## Usage

Two examples are located in `./examples`. 

* Copy the `input.txt` and the `solution.xml` the parent folder.

* Execute `flamelet`

You can change the setup in the `input.txt`.
The converged tables are stored in the folder `canteraTables`.
These can be read by the OpenFOAM utility `canteraToFoam`.

## Notes

More information on using the counterflow flame solver is available on:
https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/flameletFoam

Please feel free to contact me should you find bugs or have suggestions how to make the code better.
