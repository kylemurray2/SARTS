for updating stacks

#adjust time range, download new slcs, then adjust it back
setupStack.py

#somehow, the xml files for geom files are removed if you try to make new ones when the old ones are still in the directory.  So remove them before running isce.

rm merged/geom_reference/*.rdr
rm merged/geom_reference/*.full
rm merged/geom_reference/*.full.xml
rm merged/geom_reference/*.full.aux.xml
rm merged/geom_reference/*.full.vrt

runISCE.py

adjustGeom.py -drcf  # use the f flag if the file path has changed. 

rm -r Fringe/coreg_stack

runFringe.py -tsa # This will make new .vrt files, process the additional slcs, and adjust the stacks