#!/bin/bash

if [[ "$1" != */ ]]
then
    if [ $# -ge 1 ]
    then
	path=$1/
    else
	path=./
    fi	
else
    path=$1
fi 
rm ${path}log 2>/dev/null
rm ${path}snap* 2>/dev/null
rm ${path}fort.* 2>/dev/null
rm ${path}OUT* 2>/dev/null
rm ${path}ESC 2>/dev/null
rm ${path}raw_input 2>/dev/null
rm ${path}HIARCH 2>/dev/null
rm ${path}Bin_data.pk 2>/dev/null
rm ${path}*.hdf5 2>/dev/null
