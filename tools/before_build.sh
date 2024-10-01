#!/bin/bash

# # Check if cmake is installed
# if ! command -v cmake &> /dev/null; then
#     python3 -m pip install cmake
# else
#     echo "cmake is already installed"
# fi
#
# Install make if not already installed

if ! yum list installed make > /dev/null 2>&1; then
    yum install -y make
fi

# Install wget if not already installed
if ! yum list installed wget > /dev/null 2>&1; then
    yum install -y wget
fi

# Install OpenBLAS only if not already installed
if [ ! -d /usr/local/OpenBLAS-0.3.28 ]; then
    cd /usr/local
    wget https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.28/OpenBLAS-0.3.28.tar.gz
    tar xzf OpenBLAS-0.3.28.tar.gz
    cd OpenBLAS-0.3.28
    make DYNAMIC_ARCH=1 NUM_THREADS=64
    make install DYNAMIC_ARCH=1 NUM_THREADS=64 PREFIX=/usr/local
fi

# Install Armadillo only if not already installed
# if [ ! -d /usr/local/armadillo-14.0.2 ]; then
#     cd /usr/local
#     wget https://sourceforge.net/projects/arma/files/armadillo-14.0.2.tar.xz
#     tar xvf armadillo-14.0.2.tar.xz
#     cd armadillo-14.0.2
#     cmake . -DOPENBLAS_PROVIDES_LAPACK=true
#     make install
# fi
