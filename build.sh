#!/bin/bash

# Path to this script's directory
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CGAL_path="$scriptDir/CGAL-apps"
output_path="../../code"

echo "compiling CGAL apps"
cd $CGAL_path
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=$output_path -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE=$output_path 
make

echo "CGAL apps done!"
