#!/usr/bin/env bash
#change to the working directory
cd $1

IFS=$'\n'
for file in $(find . -name '*.gdml');
do
  mv -- "$file" "${file// /_}"
done
for file in $(find . -name '*.stl');
do
  mv -- "$file" "${file// /_}"
done
unset IFS 
