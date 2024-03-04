Wrapper for making PS/DS InSAR time series. For Linux Ubuntu  

Must first install:  
-isce2  
-Dolphin  
see docs/installNotes.sh for all commands needed to do these installs.  

Workflow:  
1. Create working directory and cd into it. Copy the params.yaml file by running:  

    get_params_sarts.py  

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

    setupStack.py

5. Run the stack processor/ISCE:

    runISCE.py

    This will likely take multiple days to run.  You can check progress by looking at logFiles/runlog_# files.   
    To monitor, you can watch the output of runlog_1 for example using:  
        watch -n 1 tail -n 30 logFiles/runlog_1  

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


9. Run Dolphin:

    At this point you could simply make downlooked IFGs with the coregistered stack and then unwrap those. Do this by running ifgs.py -dumf.  This will turn the Dolphin mode off and make the unwrapped ifgs. By default, we are doing PS-DS analysis with Dolphin, so we can run that using:

    runDolphin.py -ps  

    This runs all of the Dolphin steps to do the PS_DS analysis and should output Dolphin/adjusted_wrapped_DS/*slc  

    To rerun or run any specific step, these are the options:   

    Run sequential estimator with phase linking  

    options:  
    -h, --help           show this help message and exit  
    -p, --do_ps          Run PS estimation. (default: False)  
    -s, --sequential_PL  Run sequential_PL function. (default: False)  



10. Make interferograms, downlook, coherence, filter, and unwrap:  

    ifgs.py -dum  

    optional arguments:    
    -h, --help            show this help message and exit    
    -d, --downlook        Downlook interferograms (default: True)   
    -u, --unwrap          Unwrap interferograms (default: True)  
    -m, --make-ifgs       Make the interferograms (default: True)  
    -n, --nproc           Number of parallel processes. Use 1 for no parallelization (default: 5)   
    -f, --nodolphin        Use this flag if you are not using dolphin psds. (default: True)    


11. Recommended to use MintPy to estimate the time series and velocities. MintPy can use the unwrapped ifgs and geom_files we created. 

Install through conda/mamba:  
mamba install -c conda-forge mintpy  

It is well mantained and documented:  
https://github.com/insarlab/MintPy/tree/main  

Good example of processing workflow here:  
https://github.com/insarlab/MintPy-tutorial/blob/main/workflows/smallbaselineApp.ipynb  

To process ALOS-1:
    run get_params_sarts.py with the -a flag  
    modify params.yaml with the path, polygon, bounds, etc.  
    downloadData.py -sd -srtm (don't use orbit flag). This downloads the srtm dem   
    put the path to the .wgs84 dem in params.yaml  
    use 'prepRawALOS_modified.py -i SLCS -o SLCS' to organize/reformat the files  
    setupStack_alos.py  
    Then follow the normal steps... runISCE.py etc...  




Potential errors:  
-->    Exception: Could not determine a suitable burst offset  
fix: The orbit file may have not downloaded correctly. try redownloading it.  


notes:  
If you want to update the stack and you deleted your SLC safe files, you should adjust the search dates to be 6 acquisitions before the first date of the update so you only download the SLCs that you need.  After downloading all of those, put the search date back to the original search date before running setupStack.py (setupStack.py also uses the date range to decide what to process).    

IW mode Sentinel-1 SLCs
    Resolution: 2.7x22 m to 3.5x22 m	
    Pixel Spacing: 2.3x14.1 m

