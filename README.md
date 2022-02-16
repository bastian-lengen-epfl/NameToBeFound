# Distance ladder

This folder contains all the script to process multivariate regression on the Cepheid- and/or TRGB-based distanced ladder.

Results Blablabla

## Beforehand
Two thing have to be done before using the different functions:
* Modify the `Fit_parameters.py` in order to choose what kind of fit you want
* Include your data in the `/Data` folder. According to your need the following .csv have to be in this folder:
  * `/Data/Cepheids.csv` containing the data from Cepheids in SN-host galaxies.
  * `/Data/Cepheids_MW.csv` containing the data from the MW-Cepheids.
  * `/Data/Cepheids_anchors.csv` containing the data from the Cepheids in anchor galaxies.
  * `/Data/TRGB.csv` containing the data for Cepheids in SN-host galaxies.
  * `/Data/TRGB_anchors.csv` containing the data from the Cepheids in anchor galaxies.
  * `/Data/SNe_Cepheids.csv` containing the data from SNe in Cepheid-host galaxies.
  * `/Data/SNe_TRGB.csv` containing the data from SNe in TRGB-host galaxies.
  * `/Data/SNe_Hubble.csv` containing the data from high-z SNe for the redshift-magnitude diagram.

### Cepheids Data 
The .csv have to be in a specific form in order to run the code. For the **Cepheids** the columns have to be:
1) Galaxy host
2) log10 of the pulsation period
3) Wesenheit magnitude
4) Uncertainty on the Wesenheit magnitude
5) Metallicity [Fe/H] or [O/H] (to be decided)
6) Redshift ot the host galaxy

Additionally, for the `/Data/Cepheids_MW.csv`:
7) Parallax 
8) Uncertainty on the parallax

And for the `/Data/Cepheids_anchors.csv`:
7) Geometric distance modulus of the anchor galaxy 
8) Uncertainty on the geometric distance modulus

### TRGB Data 
For the *TRGB*, the columns have to be:
1) Galaxy host
2) Observed I-band magnitude of the TRGB
3) Uncertainty on the observed magnitude
4) Absorbtion in the I-band
5) Redshift of the host galaxy
6) A color term V-I

Additionally, for the `/Data/TRGB_anchors.csv`:
7) Geometric distance modulus of the anchor galaxy
8) Uncertainty on the geometric distance modulus

### SNe Data
For the *SN*, the columns have to be:
1) Galaxy host
2) Apparent B magnitude
3) Uncertainty on the apparent magnitude

The columns are different for the `SNe_Hubble.csv`:
1) SN name
2) Apparent B magnitude
3) Uncertainty on the apparent magnitude
4) Redshift
5) Uncertainty on the redshift

### Note
Note that utility functions that convert data from Riess 2016/2019/2021, 
Anand 2021 and the Pantheon into the desired format are available in the `/Functions/Usefull` folder

## Fit