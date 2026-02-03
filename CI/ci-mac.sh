#!/bin/bash

#Exit on any error
set -e

PATH=$PATH:/opt/local/bin
export PATH

OPTIONS=""

echo "Building with C++${CMAKE_CXX_STANDARD}";
BUILDPATH="$BUILDPATH-cpp${CMAKE_CXX_STANDARD}"

#########################################
# Build release version
#########################################

if [ ! -d build-release ]; then
  mkdir build-release
fi

cd build-release

cmake -DCMAKE_BUILD_TYPE=Release $OPTIONS ../

#build it
make

#build the unit tests
make unittests


#########################################
# Run Release Unittests
#########################################
cd Unittests

#execute tests
../Build/bin/unittests --gtest_color=yes --gtest_output=xml

cd ..
cd ..

#########################################
# Build Debug version and Unittests
#########################################

if [ ! -d build-debug ]; then
  mkdir build-debug
fi

cd build-debug

cmake -DCMAKE_BUILD_TYPE=Debug $OPTIONS ../

#build it
make

#build the unit tests
make unittests

#########################################
# Run Debug Unittests
#########################################

cd Unittests

# Run the unittests
../Build/bin/unittests --gtest_color=yes --gtest_output=xml
