#!/bin/bash
QUALIKIZ=$PWD
read -p 'Please input desired run directory name, e.g. '\''run8_JETITB'\'': ' RUNNAME
RUNPATH=$QUALIKIZ/runs/$RUNNAME
mkdir -p $QUALIKIZ/runs
mkdir $RUNPATH

ln -s $QUALIKIZ/tools/qualikiz_input.py $RUNPATH
ln -s $QUALIKIZ/QuaLiKiz $RUNPATH
