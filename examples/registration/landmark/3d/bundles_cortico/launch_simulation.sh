#!/bin/bash


#Moving to the local deformetrica/examples directory
dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $dir

# Atlas construction :
deformetrica registration 3D model.xml data_set.xml optimization_parameters.xml --output-dir=output


