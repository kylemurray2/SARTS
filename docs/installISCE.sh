#!/bin/bash

env_name="myenv"

if conda env list | grep -q "$env_name"; then
  echo "Conda environment $env_name already exists"
else
  echo "Creating conda environment $env_name"
  conda env create -f isce_fringe.yml
fi

source activate isce

if [ ! -d "$HOME/Software" ]; then
  mkdir Software
fi

cd $HOME/Software

if [ ! -d "$HOME/Software/isce" ]; then
  if [ ! -d "$HOME/Software/src" ]; then
    mkdir src
  fi

  cd ~/Software/src/
  git clone https://github.com/isce-framework/isce2.git

  # Make the buld directory
  cd ~/Software/src/isce2
  mkdir build  && cd build

  # Do the cmake stuff
  # with cuda
  #cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/Software/isce -DPYTHON_MODULE_DIR=isce2 -DCMAKE_CUDA_FLAGS="-arch=sm_75" -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} -DCMAKE_BUILD_TYPE=Release
  #without cuda
  cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/Software/isce -DPYTHON_MODULE_DIR=$HOME/Software/isce/python -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} -DCMAKE_BUILD_TYPE=Release

# Do the rest of the make:
  make -j 16 # to use multiple threads
  make install

  echo ' add the paths to your environment'

fi

# For ISCE-2
export ISCE_ROOT=/home/km/Software/isce
export ISCE_SRC_ROOT=/home/km/Software/src/isce2
# ISCE installed applications
export PATH=$PATH:$ISCE_ROOT:$ISCE_ROOT/bin:$ISCE_ROOT/applications
# ISCE applications from src for PATH
export PATH=$PATH:$ISCE_SRC_ROOT/contrib:$ISCE_SRC_ROOT/contrib/stack/topsStack:$ISCE_SRC_ROOT/contrib/stack/stripmapStack
export PATH=$PATH:$ISCE_SRC_ROOT/bin:$ISCE_SRC_ROOT/applications
# ISCE applications from src for pythonpath
export PYTHONPATH=$PYTHONPATH:$HOME/Software:$ISCE_SRC_ROOT/applications:$ISCE_SRC_ROOT/contrib:$ISCE_SRC_ROOT/contrib/stack
export PYTHONPATH=$PYTHONPATH:$ISCE_SRC_ROOT/contrib/stack/topsStack
# ISCE applications from install for pythonpath
export PYTHONPATH=$PYTHONPATH:$ISCE_ROOT
export PYTHONPATH=$PYTHONPATH:$ISCE_ROOT/components
