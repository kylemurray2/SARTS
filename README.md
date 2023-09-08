Wrapper for making PS/DS InSAR time series.

Create Conda library with docs/isce_fringe.yml file.
conda env create -f docs/isce_fringe.yml

Must first install:
-isce2
-FRINGE (forked version from kylemurray2/)

Workflow:
1. Copy the localParams_template.py to localParams.py in your working directory.

2. Edit that file with all of the settings you need, like bounding lat/lon, etc.

3. Download SLCS and Orbits and DEM:
    downloadData.py

    *Note: You need an account with Earthdata for this script to work. Make an account here:
    https://urs.earthdata.nasa.gov/
    Go to "applications" in your earth data account and approve ASF data access.
    Then save you login credentials by editing SARTS/docs/.netrc with your account info and move it to your home directory.

    *Note: To download the DEM (required), you need to get an API key from OpenTopography. Make a file in your home directory called '~/.otkey' containing just the key on the first line. The script getDEM.py is run in downloadData.py and will download the the Copernicus DEM. Copernicus DEM is recommended if processing more recent data. 
    
    Alternatively, you can get the 30 m SRTM DEM by using the --get-srtm flag. 

4. Check that the SLCs can be opened and then setup the run_files and configs:
    setupstack.py

5. Run the stack processor/ISCE:
    runISCE.py

6. Check that merged/SLC/*/*full has all of the full SLCs.  

7. Add the crop bounds to localParams.py


8. Setup Fringe.
    setupFringe.py

    This will crop the geometry files, do the watermask, and save some stuff
    to the ps namespace.  

9. Run Fringe
    runFringe.py
    This does the PS_DS analysis and should output Fringe/adjusted_wrapped_DS/*slc

10. Make interferograms, downlook, coherence, filter, and unwrap
    ifgs.py
