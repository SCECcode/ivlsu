# IVLSU  

<a href="https://github.com/sceccode/ivlsu.git"><img src="https://github.com/sceccode/ivlsu/wiki/images/ivlsu_logo.png"></a>

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![GitHub repo size](https://img.shields.io/github/repo-size/sceccode/ivlsu)
[![ivlsu-ucvm-ci Actions Status](https://github.com/SCECcode/ivlsu/workflows/ivlsu-ucvm-ci/badge.svg
)](https://github.com/SCECcode/ivlsu/actions)

SSIP Imperial Valley model provides P-wave velocities for the shallow crust in
the southern Salton Trough, from the southern Salton Sea to the USA-Mexico border.
The model combines travel-time observations from explosions recorded during the 
Salton Seismic Imaging Project (SSIP) and the IV1979 active source experiment 
with local earthquake recordings. The model extends from the surface to 8-km depth.

Persaud P., (2016) Integrated 3-D Seismotectonic-Velocity Model of the Salton Trough Region based on Seismicity and Explosive Shots from the Salton Seismic Imaging Project Final Technical Report to U.S. Geological Survey NEHRP, Award G15AP00062

## Installation

This package is intended to be installed as part of the UCVM framework,
version 25.7 or higher. 

## Library

The library ./lib/libivlsu.a may be statically linked into any
user application. Also, if your system supports dynamic linking,
you will also have a ./lib/libivlsu.so file that can be used
for dynamic linking. The header file defining the API is located
in ./include/ivlsu.h.

## Contact the authors

If you would like to contact the authors regarding this software,
please e-mail software@scec.org. Note this e-mail address should
be used for questions regarding the software itself (e.g. how
do I link the library properly?). Questions regarding the model's
science (e.g. on what paper is the IMPERIAL valley velocity model
based?) should be directed to the model's authors, located in the
AUTHORS file.


