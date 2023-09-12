Wrapper for making PS/DS InSAR time series.

Must first install:
-isce2
-FRINGE (forked version from kylemurray2/)
see docs/installNotes.sh for all commands needed to do these installs. 

Workflow:
1. Create working directory and cd into it. Copy the params.yaml file by running:

    getParams.py

2. Edit the params.yaml file with all of the settings you need. Spend some extra time double-checking all of these parameters are exactly how you want them before moving on, because most issues will be related to mistakes in this step. 

    -Bounding lat/lon, and a polygon of your area of interst for ASF search
    -Path number
    -Swath number(s)
    -reference date
    -Downlooking factors (azimuth and range)
    -Date ranges
    -Make sure to set the aux_cal and orbit path to a directory where you want to save those.  If they don't exists, they will be created, but it's best to always point to a central source which can be used in future stacks. 
    -Adjust any other parameters before starting.   

3. Download SLCS and Orbits and DEM:

    downloadData.py -sdo

    optional arguments:
    -h, --help            show this help message and exit
    -s, --search-data     Search ASF for data and output to out.csv (default: False)
    -d, --download-slc    download SLCs from ASF (default: False)
    -o, --download-orbits  download orbit files (default: False)
    -srtm, --get-srtm     Use SRTM dem instead of copernicus (default: False)

    *Note: You need an account with Earthdata for this script to work. Make an account here:
    https://urs.earthdata.nasa.gov/
    Go to "applications" in your earth data account and approve ASF data access.
    Then save you login credentials by editing SARTS/docs/.netrc with your account info and move it to your home directory.

    *Note: To download the copernicus DEM, you need to get an API key from OpenTopography. Make a file in your home directory called '~/.otkey' containing just the key on the first line. The script getDEM.py is run in downloadData.py and will download the the Copernicus DEM. Copernicus DEM is recommended if processing more recent data. 
    
    Alternatively, you can get the 30 m SRTM DEM by using the --get-srtm flag. 

4. Check that the SLCs can be opened and then setup the run_files and configs:

    setupstack.py

5. Run the stack processor/ISCE:

    runISCE.py

    This will likely take multiple days to run.  You can check progress by looking at logFiles/* files. 

6. If step 5 completed without errors, check that merged/SLC/*/*full has all of the full SLCs.  

7. Add the crop bounds to params.yaml

8. Setup files and params for the rest of the processing:

    adjustGeom.py -dcr

    This will crop the geometry files, make a watermask (optional), and save some stuff
    to the ps namespace.  

    -h, --help        show this help message and exit
    -d, --downlook    Downlook geometry files. (default: False)
    -c, --crop        Crop geometry files. (default: False)
    -r, --replace     Overwrite cropped/downlooked geometry files (default: False)
    -f, --fix-images  Fix file path in xml files (use if files were moved to different directory)
    -p, --plot-off    Turn plotting off (default: True)


9. Run Fringe:

    At this point you could simply make downlooked IFGs with the coregistered stack and then unwrap those. Do this by running ifgs.py -dumf.  This will turn the Fringe mode off and make the unwrapped ifgs. By default, we are doing PS-DS analysis with Fringe, so we can run that using:

    runFringe.py
    This does the PS_DS analysis and should output Fringe/adjusted_wrapped_DS/*slc

10. Make interferograms, downlook, coherence, filter, and unwrap:

    ifgs.py -dum

    optional arguments:
    -h, --help            show this help message and exit
    -d, --downlook        Downlook interferograms (default: True)
    -u, --unwrap          Unwrap interferograms (default: True)
    -m, --make-ifgs       Make the interferograms (default: True)
    -n NUM_PROCESSES, --nproc NUM_PROCESSES
                            Number of parallel processes. Use 1 for no parallelization (default: 5)
    -f, --noFringe        Use this flag if you are not using Fringe psds. (default: True)


11. Recommended to use MintPy to estimate the time series and velocities. MintPy can use the unwrapped ifgs and geom_files we created. 

Install through conda/mamba:
mamba install -c conda-forge mintpy

It is well mantained and documented:
https://github.com/insarlab/MintPy/tree/main

Good example of processing workflow here:
https://github.com/insarlab/MintPy-tutorial/blob/main/workflows/smallbaselineApp.ipynb

