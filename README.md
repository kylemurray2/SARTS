Wrapper for making PS/DS InSAR time series.

Create Conda library with docs/isce_fringe.yml file.

Must first install:
-isce2
-FRINGE (forked version from kylemurray2/)

Workflow:
1. Copy the localParams_template.py to localParams.py in your working directory.

2. Edit that file with all of the settings you need, like bounding lat/lon, etc.

3. Download SLCS and Orbits
    downloadData.py

4. Download SAR SLCs and Orbit files and DEM
    downloadData.py

5. Check that the SLCs can be opened and setup the run_files and configs
    runsetupStack.

6. Run the stack processor/ISCE
    runISCE.py

7. Setup Fringe.
    setupFringe.py

    This will crop the geometry files, do the watermask, and save some stuff
    to the ps namespace.  

8. Run Fringe
    runFringe.py
    This does the PS_DS analysis and should output Fringe/adjusted_wrapped_DS/*slc

9. Make interferograms, downlook, coherence, filter, and unwrap
    ifgs.py
