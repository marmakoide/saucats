#!/usr/bin/env sh

# Installs eigen 3.3.7 to $HOME/local

set -ex

wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2
tar xjf 3.3.7.tar.bz2
rm 3.3.7.tar.bz2
cd eigen-eigen-323c052e1731
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/local
make install
