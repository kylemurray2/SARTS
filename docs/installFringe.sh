#!/bin/bash

env_name="isce"

if conda env list | grep -q "$env_name"; then
  echo "Conda environment $env_name already exists"
else
  echo "Creating conda environment $env_name"
  conda env create -f $HOME/Software/SARTS/docs/isce_fringe.yml
fi

source activate isce
cd $HOME/Software

if [ ! -d "$HOME/Software/Fringe" ]; then
  mkdir Fringe
  cd Fringe
  mkdir install build src
  cd src
  git clone https://github.com/kylemurray2/fringe.git
  cd ../build
  CXX=${CXX} cmake -DCMAKE_INSTALL_PREFIX=../install ../src/fringe
  make all
  make install
  echo ' '
  echo ' Add the following to your environment:'
  echo 'export PATH=$PATH:$HOME/Software/Fringe/install/bin'
  echo 'export PYTHONPATH=$PYTHONPATH:$HOME/Software/Fringe/install/bin'
  echo 'export LD_PRELOAD=$HOME/Software/miniconda3/envs/isce/lib/libmkl_core.so.2:$HOME/Software/miniconda3/envs/isce/lib/libmkl_sequential.so.2:$HOME/Software/miniconda3/envs/isce/lib/libmkl_avx512.so.2:$HOME/Software/miniconda3/envs/isce/lib/libmkl_def.so.2'
  echo ' '

fi
