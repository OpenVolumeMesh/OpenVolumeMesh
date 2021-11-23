#!/bin/bash

COMPILER=$1
LANGUAGE=$2

# Exit script on any error
set -e 

OPTIONS=""
BUILDPATH=""

if [ "$COMPILER" == "gcc" ]; then
  echo "Building with GCC";
  BUILDPATH="gcc"
elif [ "$COMPILER" == "clang" ]; then

  OPTIONS="$OPTIONS -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang"
  echo "Building with CLANG";
  BUILDPATH="clang"  
fi  

if [ "$LANGUAGE" == "C++17" ]; then
  echo "Building with C++17";
  BUILDPATH="$BUILDPATH-cpp17"  
elif [ "$LANGUAGE" == "C++20" ]; then
  echo "Building with C++20";
  OPTIONS="$OPTIONS -DOVM_CXX_STANDARD=20"
  BUILDPATH="$BUILDPATH-cpp20"  
fi

#########################################

# Make release build folder
if [ ! -d build-release-$BUILDPATH ]; then
  mkdir build-release-$BUILDPATH
fi

cd build-release-$BUILDPATH

cmake -DCMAKE_BUILD_TYPE=Release -DSTL_DEBUG=ON $OPTIONS ../

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

if [ ! -d build-debug-$BUILDPATH ]; then
  mkdir build-debug-$BUILDPATH
fi

cd build-debug-$BUILDPATH

cmake -DCMAKE_BUILD_TYPE=Debug -DSTL_DEBUG=ON $OPTIONS ../

#build the unit tests
make unittests


#########################################
# Run Debug Unittests
#########################################

cd Unittests

# Run the unittests
../Build/bin/unittests --gtest_color=yes --gtest_output=xml


