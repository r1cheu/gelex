#!/bin/bash
PROJECT_DIR="${1:-$PWD}"
MAX_THREADS="${2:-64}"
# # Check if cmake is installed
# if ! command -v cmake &> /dev/null; then
#     python3 -m pip install cmake
# else
#     echo "cmake is already installed"
# fi

# Install make if not already installed

if ! yum list installed make > /dev/null 2>&1; then
    yum install -y make
fi

# Install wget if not already installed
if ! yum list installed wget > /dev/null 2>&1; then
    yum install -y wget
fi

# # Install OpenBLAS only if not already installed
if [[ "$BUILD_OPENBLAS" = "true" ]]; then
    echo "Building OpenBLAS"
    if [ ! -d /usr/local/OpenBLAS-0.3.28 ]; then
        cd /usr/local
        wget https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.28/OpenBLAS-0.3.28.tar.gz
        tar xzf OpenBLAS-0.3.28.tar.gz
        cd OpenBLAS-0.3.28
        make DYNAMIC_ARCH=1 NUM_THREADS=${MAX_THREADS} USE_OPENMP=1 QUIET_MAKE=1 NO_AFFINITY=1 USE64BITINT=1 TARGET=HASWELL,SKYLAKEX USE_LAPCAK=1 CROSS=1
        make install DYNAMIC_ARCH=1 NUM_THREADS=${MAX_THREADS} USE_OPENMP=1 QUIET_MAKE=1 NO_AFFINITY=1 USE64BITINT=1 TARGET=HASWELL,SKYLAKEX CROSS=1 USE_LAPACK=1 PREFIX=/usr/local
    fi
fi

# if [[ "$INSTALL_OPENBLAS" = "true" ]] ; then
#     echo "Installing OpenBLAS"
#     echo PKG_CONFIG_PATH $PKG_CONFIG_PATH
#     PKG_CONFIG_PATH=$PROJECT_DIR/.openblas
#     echo "${PKG_CONFIG_PATH}"
#     echo "${PROJECT_DIR}"
#     rm -rf $PKG_CONFIG_PATH
#     mkdir -p $PKG_CONFIG_PATH
#     python -m pip install scipy-openblas64
#     export LD_LIBRARY_PATH=$PKG_CONFIG_PATH/lib:$LD_LIBRARY_PATH
#     # Copy the OpenBLAS libraries to the PKG_CONFIG_PATH and rename the libscipy_openblas64.so to libopenblas.so
#     python <<EOF
# import os, scipy_openblas64, shutil
# srcdir = os.path.join(os.path.dirname(scipy_openblas64.__file__), "lib")
# shutil.copytree(srcdir, os.path.join("$PKG_CONFIG_PATH", "lib"))
# EOF
#     # pkg-config scipy-openblas --print-provides
# fi
#
# Install Armadillo only if not already installed
# if [ ! -d /usr/local/armadillo-14.0.2 ]; then
#     cd /usr/local
#     wget https://sourceforge.net/projects/arma/files/armadillo-14.0.2.tar.xz
#     tar xvf armadillo-14.0.2.tar.xz
#     cd armadillo-14.0.2
#     cmake . -DOPENBLAS_PROVIDES_LAPACK=true
#     make install
# fi
