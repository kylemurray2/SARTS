#!/bin/bash

# This file contains the commands needed to install isce, fringe, mintpy, and sarts. 
#   It also contains the environment setup commands

#Define path where you want the software to be located
softwareDir=$HOME/Software
cd $softwareDir

# Git SARTS. This has a requirements file we'll use 
# git clone https://github.com/kylemurray2/SARTS.git

# First get mamba if it doesn't exist
if [ -d $softwareDir/mambaforge ]; then
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
    chmod +x Mambaforge-Linux-x86_64.sh
    ./Mambaforge-Linux-x86_64.sh -b -p $softwareDir/mambaforge
    rm Mambaforge-Linux-x86_64.sh
fi

export PATH=$softwareDir/mambaforge/bin:$PATH

# Make the conda env that will work for isce, fringe, SARTS, and mintpy. Takes awhile..
mamba update -n base -c conda-forge mamba -y
mamba env create -f $softwareDir/SARTS/docs/requirements.yml
source activate isce

# Install ISCE-2
cd $softwareDir
mkdir -p src
if [ -f $softwareDir/src/isce2 ]; then
    rm -rf $softwareDir/src/isce2
fi
cd $softwareDir/src/
git clone https://github.com/isce-framework/isce2.git
cd $softwareDir/src/isce2
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$softwareDir/isce2 -DPYTHON_MODULE_DIR=$softwareDir -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} -DCMAKE_BUILD_TYPE=Release
make -j 10 # to use multiple threads
make install


# Install Fringe
cd $softwareDir
if [ -f $softwareDir/Fringe ]; then
    rm -rf $softwareDir/Fringe
fi
mkdir -p $softwareDir/Fringe
cd $softwareDir/Fringe
mkdir install build src
cd src
git clone https://github.com/kylemurray2/fringe.git
cd ../build
CXX=${CXX} cmake -DCMAKE_INSTALL_PREFIX=../install ../src/fringe
make all
make install

# Add the following to a setup file and source that file to export the environments
#__________________________________________________________________________________
#----------------------------------------------------------------------------------
softwareDir=$HOME/Software
# mambaforge
export PATH=$PATH:$softwareDir/mambaforge/bin

# Fringe
export PATH=$PATH:$softwareDir/Fringe/install/bin
#export LD_PRELOAD=$softwareDir/mambaforge/envs/isce/lib/libmkl_core.so:$softwareDir/mambaforge/envs/isce/lib/libmkl_sequential.so:$softwareDir/mambaforge/envs/isce/lib/libmkl_avx512.so:$softwareDir/mambaforge/envs/isce/lib/libmkl_def.so
export PYTHONPATH=$PYTHONPATH:$softwareDir/Fringe/install/bin
export PYTHONPATH=$PYTHONPATH:$softwareDir/Fringe/install/python

# ISCE
export ISCE_ROOT=$softwareDir/isce
export ISCE_SRC_ROOT=$softwareDir/src/isce2

export PATH=$PATH:$ISCE_ROOT:$ISCE_ROOT/bin:$ISCE_ROOT/applications
export PATH=$PATH:$ISCE_SRC_ROOT/contrib/stack/topsStack

export PYTHONPATH=$PYTHONPATH:$ISCE_ROOT
export PYTHONPATH=$PYTHONPATH:$ISCE_ROOT/components
export PYTHONPATH=$PYTHONPATH:$softwareDir:$ISCE_SRC_ROOT/applications:$ISCE_SRC_ROOT/contrib:$ISCE_SRC_ROOT/contrib/stack
export PYTHONPATH=$PYTHONPATH:$ISCE_SRC_ROOT/contrib/stack/topsStack

# SARTS
export PATH=$PATH:$softwareDir/SARTS
export PYTHONPATH=$PYTHONPATH:$softwareDir/SARTS

source activate isce

#----------------------------------------------------------------------------------
#__________________________________________________________________________________

# Potential errors:

#   ImportError: libgdal.so.31: cannot open shared object file: No such file or directory
# Solution:
#   ln -s $softwareDir/mambaforge/envs/isce/lib/libgdal.so.33 $softwareDir/mambaforge/envs/isce/lib/libgdal.so.31

