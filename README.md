## SpekQual

**SpekQual** is a Python script for X-ray spectrum quality assessment, with this
code one can determine beam quality characteristics, such as: HVL<sub>1</sub>, HVL<sub>2</sub>,
QVL, TVL, H<sub>i</sub>, etc. Currently spectra data from [FLUKA](https://fluka.cern/)
output and [SpekPy](https://bitbucket.org/spekpy/spekpy_release/wiki/Home)
are supported, other formats could be added in the future releases.

### Quick installation guide

Install latest Python 3 interpreter and pip tool, then clone the repository to the
prefered location:

    git clone https://github.com/GordoNice/spek_qual.git

This command will download all the files in current directory.

### Requirements

Additional modules for SpekQua; use (automatically install via pip tool):

* numpy>=1.21.0
* pandas>=1.5.3
* scipy>=1.3.3

All modules listed in `requirements.txt`.

### Basic usage

Script works from cli, type:

    ./spek_qual.py -h

to show help:

    usage: spek_qual.py [-h] [-f format] [-fdn fluka_det_n] [-frd fluka_row_drop] [-fs] [-ss spekpy_sep] [-ms mu_source] [-sc] [-v] spec.dat
    
    spek_qual script calculates beam quality characteristics for the specified X-ray spectrum
    
    positional arguments:
      spec.dat              Input spectrum file
    
    optional arguments:
      -h, --help            show this help message and exit
      -f format, --format format
                            Format of the spectrum file data: "fluka", "spekpy" (default: fluka)
      -fdn fluka_det_n, --fluka_det_n fluka_det_n
                            Only for the spectrum data format: detector number in file to load (Detector n) (default: 1)
      -frd fluka_row_drop, --fluka_row_drop fluka_row_drop
                            Only for the FLUKA spectrum data format: how many rows to drop from beginning of the file (cut to 1 keV) (default: 0)
      -fs, --fluka_save     Only for the FLUKA spectrum data format: save the extracted spectrum data from FLUKA file to the "result" foolder (default: False)
      -ss spekpy_sep, --spekpy_sep spekpy_sep
                            Only for the SpekPy spectrum data format: delimiter used to separate columns in file with spectrum (default: )
      -ms mu_source, --mu_source mu_source
                            Choose mu data: "nist" or "pene" (default: nist)
      -sc, --spekpy_char    Only for the SpekPy spectrum data format: flag to get column with the characteristic fluence, otherwise total fluence will be used for the beam qualification (default: False)
      -v, --verbose         Print information while execution (default: False)

Some examples of spectra data are available under `./data/spectra` directory, all
the spectra must be stored under this path! For instance, we can load FLUKA data:

    ./spek_qual.py -f "fluka" -fdn 7 -frd 2 60keV_20_W_Al0-5_21_tab.lis

and get the results as:

    spek_qual total beam quality characteristics for the spectrum 60keV_20_W_Al0-5_21_tab.lis

    HVL1 (Al / Cu): 0.758746 / 0.0241668 mm
    HVL2 (Al / Cu): 1.3664 / 0.0432576 mm
    QVL (Al / Cu):  2.12515 / 0.0674244 mm
    TVL (Al / Cu):  5.08868 / 0.167239 mm
    hi (Al / Cu):   0.555288 / 0.558671 mm
    Eeff (Al / Cu): 20.1202 / 20.3921 keV
    Emean:  30.0545 keV

In the file `60keV_20_W_Al0-5_21_tab.lis` we have 7 spectra of differential yield
on every stage of X-ray  passage through material. Spectrum is loaded for the final stage (`-fdn 7`).
We also need to be sure that no data with energy less than 1 keV will be loaded, that's why
we need to drop the first 2 rows from the beginning of data (`-frd 2`).

All the results will be saved under `./result/60keV_20_W_Al0-5_21_tab_spekqual.dat` path.

SpekPy spectrum could be loaded as:

    ./spek_qual.py -f "spekpy" 130keV_20_W_Al0-5.dat

results will be:

    spek_qual total beam quality characteristics for the spectrum 130keV_20_W_Al0-5.dat

    HVL1 (Al / Cu): 1.47345 / 0.0490064 mm
    HVL2 (Al / Cu): 3.69617 / 0.151569 mm
    QVL (Al / Cu):  5.16962 / 0.200575 mm
    TVL (Al / Cu):  14.1183 / 0.75245 mm
    hi (Al / Cu):   0.398642 / 0.323328 mm
    Eeff (Al / Cu): 25.6106 / 26.2829 keV
    Emean:  48.8191 keV

NIST data for all the mu coefficients will be used by default! Penelope data could
be used with: `-mu_source "pene"`.

### Documentation

Will be available soon...

### TODO list
- [ ] Detailed documentation
- [ ] Nice pdf output (using LaTeX) with all the results for loaded spectrum
- [ ] Other spectra data formats (please contact me for any suggestions)
- [ ] ???
