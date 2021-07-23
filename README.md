These scripts were used for CoMetChem data analysis of the peptide H3(18-26).

# Requirements

These scripts depend on the following packages, which can be installe with `pip`
prior to running any of the scripts in this repository:

- pandas
- pyteomics
- matplotlib
- numpy
- scipy
- wx
- pastaq


# Usage
Fragmentation.py was used to generate a inclusion list for PASTAQ.
MS2_abundance_H3.py and Site-specific_intensities_H3.py were executed consecutively. The resulting output was used for further analysis.

Before running these scripts some things have to be specified:

## Fragmentation.py
- sequence of your peptide (one-letter code)
- the possible modifications (DICT_MODIFICATIONS)
- which aminoacid should be modified (line 272)

## MS2_abundance_H3.py
- path where the folders containing the data are stored (PASTAQ output)
- Dictionary converting the date to a replicate number
- Dictionary converting the timepoints contained in the folder name to numbers

## Site-specific_intensities_H3.py
- Path of the corrected MS1 data (*.csv)
