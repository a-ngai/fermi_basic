# Quick start-up (for Windows)

1. Create a python virtual environment

> python -m venv \<name\>

2. Activate the virtual environment (on Windows; Linux has a different file structure in the virtual environment)

> \<name\>\Scripts\activate

3. Within the virtual environment, install the following modules:

> pip install h5py lmfit matplotlib numpy scipy

4. Run the script

> python electron_magnetic_bottle.py

# Quick explanation

The "electron_magnetic_bottle.py" script is a simple script which performs the
following:

1. extracts magnetic bottle data from an .h5 file,
2. separates foreground from background data, based on the
 "BACKGROUND_PERIOD"  found in the .h5 file, and
3. gives a crude calibration from TOF to eKE coordinates, based on the TOF
 and eKE values given.
4. rebins spectra (you will use this regularly!).
    
This gives a simple format from which more things can be added, and additional
data from the .h5 files can be loaded.

# Python modules used in testing:

- h5py==3.12.1
- lmfit==1.3.2
- matplotlib==3.10.0
- numpy==2.2.2
- scipy==1.15.1

# Difference from "https://github.com/a-ngai/SolvatedElectronsBeamtime"

This small script is the minimum needed to read in and manipulate data saved from the detectors at FERMI. The additional functionality of the "SolvatedElectronsBeamtime" repository is as follows:

- automating extraction of fore-/background data (based on the background-period of the gas nozzle, and possibly the SLU running not at 50Hz),
- filtering data using search line operators,
- defining spectrometer calibration functions in the general case,
- quick caching and retrieval of processed data.

Although this repository has gone through testing, it is not very stable to missing data or file structure changes. It can tolerate some of these deficiencies, but thorough testing of its stability has not been done. The code structure is also unfortunately complicated, so debugging those problems may require knowledge of the repository, which may prove challenging.