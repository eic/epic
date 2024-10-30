#!/usr/bin/env bash

# Set the field separator to newline to handle filenames with spaces
IFS=$'\n'

#change to the diretory in which files are processed
cd $1
mkdir ASCII_STL

# Process .stl files
for f in $(find . -type f -name "*.stl" ! -name "*_ASCII.stl"); do
    # Unset IFS to avoid issues within the loop
    unset IFS
    # Check if the STL file is binary or ASCII
    if grep -q 'solid' "$f"; then
        echo "$f is already an ASCII STL and should be with name ending _ASCII.stl ."
        # Process the STL file with stl_gdml.py
    else
        echo "Converting $f from binary to ASCII STL."
        /cvmfs/sft.cern.ch/lcg/releases/Python/2.7.9.p1-df007/x86_64-slc6-gcc49-opt/bin/python2 ../BinaryToASCII.py "$f"    
    fi
    # Reset IFS for the next iteration
    IFS=$'\n'
done

# Process the ASCII .stl files
for f in $(find . -name '*_ASCII.stl'); do
    # Unset IFS to avoid issues within the loop
    unset IFS
    # Check if the STL file is binary or ASCII
    if grep -q 'solid' "$f"; then
        
        echo "Converting $f from ASCII STL to .gdml"
        python ../stl_gdml.py ../out.gdml "$f"
    else
        echo "$f is probably a binary STL without ASCII start token solid and should not have the _ASCII.stl name ending."  
    fi
    # Reset IFS for the next iteration
    IFS=$'\n'
done

# Unset IFS and disable filename expansion
unset IFS
set +f
