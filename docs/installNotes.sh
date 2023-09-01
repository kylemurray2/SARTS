miniconda
Fringe
MintPy
SARTS
ISCE

git clone https://github.com/kylemurray2/fringe.git
git clone https://github.com/isce-framework/isce2.git

# First get miniconda
cd $HOME/Software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
rm Miniconda3-latest-Linux-x86_64.sh

# Git SARTS. This has a requirements file we'll use 
git clone https://github.com/kylemurray2/SARTS.git
# Make the conda env that will work for isce, fringe, SARTS, and mintpy. Takes awhile..go get a snack
conda env create -f SARTS/docs/isce_fringe.yml

# Install ISCE-2
mkdir src
cd $HOME/Software/src/
git clone https://github.com/isce-framework/isce2.git
cd $HOME/Software/src/isce2
mkdir build
cd build

cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/Software/isce -DPYTHON_MODULE_DIR=$HOME/Software/isce/python -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} -DCMAKE_BUILD_TYPE=Release
make -j 16 # to use multiple threads
make install