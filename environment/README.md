# Software environment

The portable public scripts were verified with R 4.4.3 and the R packages listed in `R_packages.txt`, and with Python 3.11.15 and the packages pinned in `requirements.txt`.

The historical CPTAC analysis used Python 3.12.13, NumPy 2.3.5 and pandas 3.0.1. Historical SciPy and Matplotlib versions were not recorded and are therefore not asserted here. Re-execution in the current verified environment reproduced the retained CPTAC results within the predefined floating-point tolerance (absolute `1e-10`, relative `1e-8`).

