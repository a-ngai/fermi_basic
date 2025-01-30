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

Python modules used in testing:

h5py==3.12.1
lmfit==1.3.2
matplotlib==3.10.0
numpy==2.2.2
scipy==1.15.1