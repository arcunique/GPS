[![Build Status](https://img.shields.io/badge/release-1.0.0-orange)](https://github.com/arcunique/Cplotter)
[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-371/)

GPS (Genesis Population Synthesis)  
======

This is a Python codebase for developing population synthesis models by using the Genesis database of planet formation.
Genesis database is a database of planet formation models for small exoplanets (super-Earths and Mini-Neptunes) developed by 
[Mulders et al. 2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...897...72M/abstract) by running N-body simulations. This is a part of the 
[Alien Earths](https://eos-nexus.org/genesis-database/) project where further details about the models can be found. The data is provided with the code. However, if lost, the data can be found at [Github](https://github.com/GijsMulders/Genesis).

[GPS](https://github.com/arcunique/GPS) presently works as a suite of post-processing codes that can independently be applied
to Genesis models and produce models for the final architectures of the small exoplanetary systems. It includes codes to compute the
bulk compositions of the planets and to simulate atmospheric loss and evolution to find out the final states of the planets that can be
observationally verified. The code also offers tools to process and analyze the data from recentmost observations of small exoplanets, in 
order to compare them with the models. This code has been used to generate results for our paper [Chakrabarty & Mulders 2023]() which
focuses on the enigmatic world of *water planets* and their possible locations. The results shown in the paper can be reproduced by running the Python 
scripts placed in the *chakrabarty_mulders_2023* folder. Please refer to the release *1.0.0* to access the version that matches the results.

Follow this for more updates. If you find this code helpful in your work, please cite the following papers:
[Chakrabarty & Mulders 2023]() and [Mulders et al. 2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...897...72M/abstract). The citation details can be found on those links.


Author
------
* Aritra Chakrabarty (Data Observatory & Universidad Adolfo Ibanez, Santiago, Chile)

Requirements
------------
* python>3.6
* numpy
* pandas
* matplotlib 
* scipy
* tqdm (optional)

Instructions on installation and use
------------------------------------
Presently, the code is only available on [Github](https://github.com/arcunique/GPS). Download the zip folder, place anywhere, and
start using.

The code contains classes and functions to perform the different stages of the models. A documentation is underway. Notebooks can be found in the 
*example_notebooks* folder which illustrate how to use these classes and functions for different use cases. Although the codebase focuses on the Genesis models, 
any other models can easily be integrated with this code. Just ensure that structure of the input file  matches with that of the Genesis file.
The structure is explained in the example notebooks. 







