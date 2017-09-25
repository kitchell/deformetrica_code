#!/bin/bash

. ../find_deformetrica.sh

# Surface matching in 2D :
deformetrica matching 2D model.xml data_set.xml optimization_parameters.xml

