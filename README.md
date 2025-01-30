The "electron_magnetic_bottle.py" script is a simple script which performs the
following:

1. extracts magnetic bottle data from an .h5 file,
2. separates foreground from background data, based on the
 "BACKGROUND_PERIOD"  found in the .h5 file, and
3. gives a crude calibration from TOF to eKE coordinates, based on the TOF
 and eKE values given.
    
The script has not been tested! This gives a simple format from which more
things can be added, and additional data from the .h5 files can be loaded.