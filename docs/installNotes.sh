# This file contains the commands needed to install isce, fringe, mintpy, and sarts. 
#   It also contains the environment setup commands

#Define path where you want the software to be located
softwareDir=$HOME/Software

# First get mamba
cd $softwareDir
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
rm Mambaforge-Linux-x86_64.sh

# Git SARTS. This has a requirements file we'll use 
git clone https://github.com/kylemurray2/SARTS.git
# Make the conda env that will work for isce, fringe, SARTS, and mintpy. Takes awhile..go get a snack
mamba update -n base -c conda-forge mamba
#conda env create -f $softwareDir/SARTS/docs/requirements.yml
# use the cloud version if you don't want mdx and spyder
mamba env create -f $softwareDir/SARTS/docs/requirements_cloud.yml

# Install ISCE-2
mkdir src
cd $softwareDir/src/
git clone https://github.com/isce-framework/isce2.git
cd $softwareDir/src/isce2
mkdir build
cd build

cmake .. -DCMAKE_INSTALL_PREFIX=$softwareDir/isce2 -DPYTHON_MODULE_DIR=$softwareDir -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} -DCMAKE_BUILD_TYPE=Release
make -j 16 # to use multiple threads
make install

#copy the contrib directory from src to the install dir, or just source to src/isce2
# cp -r $softwareDir/src/isce2/contrib $softwareDir/isce/


# Install Fringe
mkdir $softwareDir/Fringe
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
# Miniconda
export PATH=$PATH:$softwareDir/mambaforge/bin

# Fringe
export PATH=$PATH:$softwareDir/Fringe/install/bin
export LD_PRELOAD=$softwareDir/mambaforge/envs/isce/lib/libmkl_core.so:$softwareDir/mambaforge/envs/isce/lib/libmkl_sequential.so:$softwareDir/mambaforge/envs/isce/lib/libmkl_avx512.so:$softwareDir/mambaforge/envs/isce/lib/libmkl_def.so
export PYTHONPATH=$PYTHONPATH:$softwareDir/Fringe/install/bin
export PYTHONPATH=$PYTHONPATH:$softwareDir/Fringe/install/python

# ISCE
export ISCE_ROOT=$softwareDir/isce2
export ISCE_SRC_ROOT=$softwareDir/src/isce2
export PATH=$PATH:$ISCE_ROOT:$ISCE_ROOT/bin:$ISCE_ROOT/applications
export PYTHONPATH=$PYTHONPATH:$ISCE_ROOT
export PYTHONPATH=$PYTHONPATH:$ISCE_ROOT/components
export PYTHONPATH=$PYTHONPATH:$softwareDir:$ISCE_SRC_ROOT/applications:$ISCE_SRC_ROOT/contrib:$ISCE_SRC_ROOT/contrib/stack
export PYTHONPATH=$PYTHONPATH:$ISCE_SRC_ROOT/contrib/stack/topsStack

# SARTS
export PATH=$PATH:$softwareDir/SARTS

source activate isce

#----------------------------------------------------------------------------------
#__________________________________________________________________________________