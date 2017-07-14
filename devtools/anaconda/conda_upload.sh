#!/bin/bash

# Only need to change these two variables
PKG_NAME=oeommtools
USER=nividic

OS=$TRAVIS_OS_NAME-64
mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
export VERSION=`date +%Y.%m.%d`
conda build .
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l main_nigthly $CONDA_BLD_PATH/$OS/$PKG_NAME-`date +%Y.%m.%d`*.tar.bz2 --force
anaconda logout