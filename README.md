# Parity-Odd-4PCF-regions
Code for reproduction and region-based analysis of BOSS CMASS parity-odd-four-point function in Krolewski, May, Smith, Hopkins (2024)

We provide the 4PCF output text files, and the python scripts and shell scripts used to process the BOSS LSS catalogs into parity-odd 4PCF measurements using [ENCORE](https://github.com/oliverphilcox/encore). We provide shell scripts in the "encore" and "encore18" directories; these refer to the scripts used for the 10-bin runs (encore) and the 18-bin runs (encore18) respectively. These directories originally had encore installed in them; we omit the encore files as they can be found at the ENCORE github page, installed with the proper compilation flags and modifying NBIN in encore.cpp to 10 or 18 depending on the desired case.

Python scripts to process the LSS catalogs are in "data." These process both the [data](https://data.sdss.org/sas/dr12/boss/lss/) and the patchy mocks for [CMASSLOWZTOT](https://data.sdss.org/sas/dr12/boss/lss/dr12_multidark_patchy_mocks/) and [CMASS](https://www.ub.edu/bispectrum/page11.html). Scripts processing the data use a "zcut" file which is just the publicly released LSS catalog with a redshift cutoff to the desired range, 0.43 < z < 0.70.

Parity-odd four point functions for the data and mocks are in "out", organized into 10-bin (encore) and 18-bin (encore18) directories. We also distribute files of chi2 summarizing the results.

The analytic covariances needed to compute chi2 are in "cov." The computation of chi2 follows the notebooks [in the Parity-Odd-4PCF repo](https://github.com/oliverphilcox/Parity-Odd-4PCF).

Many thanks to the original authors of the parity-odd four point function papers (link): Jiamin Hou, Zachary Slepian, Robert Cahn, and Oliver Philcox. This repository is intended to allow the community to reproduce our (and their) work.
