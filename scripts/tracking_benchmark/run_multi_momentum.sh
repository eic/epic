#!/bin/bash
# Script; shyam kumar; INFN Bari, Italy
# Modified: Tuna Tasali. University of Oxford
# shyam.kumar@ba.infn.it; shyam055119@gmail.com
# first argument is the name of the folder in which the results will go

#CONFIG_VARIABLES
#                    0       1      2    3      4      5    6      7    8   9     10  11  12    13   14 
mom_array=(0.50 0.75 1.0 1.25 1.75 2.0 2.50 3.0 4.0 5.0 7.0 8.5 10.0 12.5 15.0)
#note that these variables must match between the single and multi momentum scripts
particle_array=("pi-")
filename=("tracking_output") 
etabin_array=(-3.5 -2.5 -1.0 1.0 2.5 3.5)
nevents=10000
results_dir="results_test"
compact_file_name="epic_craterlake_tracking_only.xml"

if [ -d "$results_dir" ]; then
    echo "Directory '$results_dir' already exists. Aborting program to avoid overwriting data."
    exit 1
else
    mkdir $results_dir
    echo "Directory '$results_dir' is created."
fi 
mkdir -p $results_dir/truthseed/pi-/mom $results_dir/realseed/pi-/mom $results_dir/truthseed/pi-/dca $results_dir/realseed/pi-/dca

#create a new array which has the momenta formatted to two decimal places
decimal_places=2
formatted_mom_array=()
for value in "${mom_array[@]}"; do
    formatted_mom=$(printf "%.${decimal_places}f" "$value")
    formatted_mom_array+=("$formatted_mom")
done

# run the simulation and the reconstruction

for ((i=0; i<${#formatted_mom_array[@]}; i++)); do
    echo "Submit run $i '${formatted_mom_array[i]}' to cluster"
    condor_submit MOM=${formatted_mom_array[i]} run_simulation.submit
done
