#!/bin/sh

BASE_DIR="\/home\/svanhulle\/VIV\/"
# BASE_DIR="\/home\/simon\/Documents\/Delft\/HPB-Delft\/HPB-VIV\/openFoam-VIV"

dirname="U9"
magUinf=1.8
viscosity=0.018
purgewrite=3

rex_spaces="[ \t]*?"
rex_float="[0-9]+\.?[0-9]*"
tab="    "


# Velocity
perl -i -0pe "s/internalField$rex_spaces uniform$rex_spaces\($rex_float/internalField uniform ($magUinf/g" 0/U

perl -i -0pe "s/type$rex_spaces fixedValue;\n$rex_spaces value$rex_spaces uniform$rex_spaces\($rex_float/type fixedValue;\n${tab}${tab}value${tab}uniform ($magUinf/g" 0/U


# Submission script
perl -i -0pe "s/\#PBS -N [\w-]+/\#PBS -N ${dirname}/g" sub-serial.pbs

perl -i -0pe "s/cd .*/cd $BASE_DIR\/$dirname/g" sub-serial.pbs


# Viscosity
perl -i -0pe "s/nu${rex_spaces}nu${rex_spaces}\[ 0 2 -1 0 0 0 0 \]${rex_spaces}$rex_float/nu${tab}${tab}nu \[ 0 2 -1 0 0 0 0 \] $viscosity/g" constant/transportProperties


# ControlDict
perl -i -0pe "s/magUInf${rex_spaces}${rex_float}/magUInf $magUinf/g" system/controlDict

perl -i -0pe "s/purgeWrite${rex_spaces}${rex_float}/purgeWrite $purgewrite/g" system/controlDict
