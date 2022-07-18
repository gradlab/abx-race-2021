# Racial/ethnic disparities in antibiotic use and healthcare utilization

[![DOI](https://zenodo.org/badge/432842870.svg)](https://zenodo.org/badge/latestdoi/432842870)

This repository includes code needed to reproduce the results from the manuscript "Racial/ethnic disparities in antibiotic use and healthcare utilization, United States, 2016/2018: a cross-sectional study", available on [medRxiv](https://www.medrxiv.org/content/10.1101/2021.12.09.21266965v3).

## File structure

- `download-data.R`: Downloads Stata-formatted data from NAMCS/NHAMCS website; downloads ICD-10 codes from Chua et al.; saves data as `.rds` format
- `data/`: Folder for all the saved data
- `analyze.R`: Analyze downloaded data

## Authors

Scott W. Olesen <olesen@hsph.harvard.edu> <br>
Stephen M. Kissler <skissler@hsph.harvard.edu>

