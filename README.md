effective Surface Quasi-Geostrophic (eSQG)
================================================

This software applies QG framework to retrieve ocean interior state using Sea Surface Height (SSH) anomaly fields and ocean stratification as the input.

Installation
-------------

To install eSQG, simply clone this repository and run the setup file using the following command:
    
    $ python setup.py install

Dependencies
-------------

eSQG requires the following dependencies to be installed:

"numpy", "seawater", "sklearn"

Usage
-------------

eSQG can be used to reconstruct the three-dimensional subsurface ocean potential density and velocity based on SSH anomaly, effective ocean stratification climatology, and a rms amplitude of the vertical velocity.

Contributing
-------------

Contributions to eSQG are welcome! Please feel free to fork this repository and submit pull requests with your improvements.

Citation
-------------

This software is based on the effective Surface Quasi-geostrophic (isQG) framework proposed by Qiu et al. in 2016. 
