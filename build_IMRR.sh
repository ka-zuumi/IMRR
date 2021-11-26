#!/bin/bash

FC=gfortran

###########################################################################

cd bin

$FC -free -c DOUBLE_v.f90
$FC -free -c PARAMETERS_v.f90
$FC -free -c VARIABLES_v.f90
$FC -free -c FUNCTIONS_v.f90
$FC -free -c ANALYSIS_v.f90
$FC -free -c ls_rmsd_original_v.f90
$FC -free -c SIMILARITY_v.f90
$FC -free -c interactMultipleGrids_v.f90


$FC -free DOUBLE_v.o PARAMETERS_v.o VARIABLES_v.o FUNCTIONS_v.o ANALYSIS_v.o ls_rmsd_original_v.o SIMILARITY_v.o interactMultipleGrids_v.o example1.f90 -o example1.o

cd ..
