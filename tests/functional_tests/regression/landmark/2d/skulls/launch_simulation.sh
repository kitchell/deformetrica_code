#!/bin/bash

. ../find_deformetrica.sh

# Regression construction :
deformetrica regression 2D model.xml data_set.xml optimization_parameters.xml
