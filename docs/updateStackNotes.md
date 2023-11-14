# Notes for updating stacks

# adjust time range in params.yaml from the sixth most recent acquisition to the new end date and download new slcs
downloadData.py -sdo

# Then you need to adjust the start time in params.yaml back to the original, before running setupStack.py. Also remove the run_files directory
rm -r run_files
setupStack.py

# somehow, the xml files for geom files are removed if you try to make new ones when the old ones are still in the directory.  So just remove geom_reference, and we'll make new files. This could be fixed later..
rm -r merged/geom_reference

runISCE.py

adjustGeom.py -drc -f  # use the f flag if the file path has changed. 

# remove the fringe stack so we can make an updated one
rm -r Fringe/coreg_stack

# run fringe to process new SLCs. make sure the stack size is the same as it was for the original run.
runFringe.py -tsa # This will make new .vrt files, process the additional slcs, and adjust the stacks



reference
Fringe/adj
merged/baseline
merged/geom_reference
