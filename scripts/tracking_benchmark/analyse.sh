#!/bin/bash
# use: ./analyse.sh <results_dir>

# source epic/install/bin/thisepic.sh
#                    0       1      2    3      4      5    6      7    8   9     10  11  12    13   14 
mom_array=(0.50 0.75 1.0 1.25 1.75 2.0 2.50 3.0 4.0 5.0 7.0 8.5 10.0 12.5 15.0)
particle_array=("pi-")
filename=("tracking_output") 
etabin_array=(-3.5 -2.5 -1.0 1.0 2.5 3.5)
nevents=10000

#the analysis directory is the first argument
results_dir=$1
echo "Analayse files in $1"

#create a new array which has the momenta formatted to two decimal places
#also arrange the entries into an array with curly braces as that is what the root macro likes
decimal_places=2
formatted_mom_array='{'
for value in "${mom_array[@]}"; do
    formatted_mom=$(printf "%.${decimal_places}f" "$value")
    formatted_mom_array+="$formatted_mom"
    formatted_mom_array+=','
done
formatted_mom_array+='}'
echo ${formatted_mom_array[@]}

#add analysis files if they are not already in $results_dir.
analysis_files=("Tracking_Performances.C" "doCompare_truth_real_widebins_mom.C" "doCompare_truth_real_widebins_dcaz.C" "doCompare_truth_real_widebins_dcaT.C")
for file in "${analysis_files[@]}"; do
    if [ -e "$results_dir/$file" ]; then
        echo "INFO: $file already exists, not altering for analysis integrity purposes."
    else
        echo "INFO: $file does not exist in $results_dir. Copying..."
        cp "$file" "$results_dir/"
    fi
done

#change to the particular results dir to analyse the data
cd $results_dir
rm -r truthseed
rm -r realseed
mkdir -p truthseed/pi-/mom realseed/pi-/mom truthseed/pi-/dca realseed/pi-/dca

# run the analysis
for ((iparticle=0; iparticle<${#particle_array[@]}; iparticle++)); do
# truth seeding
for ((i=0; i<${#mom_array[@]}; i++)); do
#Form("./%s_%1.1f",mom)
root -b -l -q Tracking_Performances.C'("'${filename}_$(printf "%.2f" ${mom_array[i]})'.edm4eic.root","'${particle_array[iparticle]}'",'${mom_array[i]}',0.15,false)'
done

# real seeding
for ((i=0; i<${#mom_array[@]}; i++)); do
root -b -l -q Tracking_Performances.C'("'${filename}_$(printf "%.2f" ${mom_array[i]})'.edm4eic.root","'${particle_array[iparticle]}'",'${mom_array[i]}',0.15,true)'
done
done
cd truthseed/pi-/dca
hadd final_hist_dca_truthseed.root *.root
cd ../../../realseed/pi-/dca/
hadd final_hist_dca_realseed.root *.root
cd ../../../


rm -rf Final_Results/ Debug_Plots/ 
mkdir -p Final_Results/pi-/mom Final_Results/pi-/dca  Debug_Plots/truth/pi-/mom  Debug_Plots/truth/pi-/dca  Debug_Plots/real/pi-/mom Debug_Plots/real/pi-/dca
# loop over particles
for ((iparticle=0; iparticle<${#particle_array[@]}; iparticle++)); do
for ((i=0; i<${#etabin_array[@]}-1; i++)); do
xmax_hist=0.3 
if [ $i == 2 ] || [ $i == 1 ]; then 
xmax_hist=0.01 
fi

#print the momenta into an array with curly braces

root -b -l -q doCompare_truth_real_widebins_mom.C'("'${particle_array[iparticle]}'",'${etabin_array[i]}','${etabin_array[i+1]}','$xmax_hist')'
root -b -l -q doCompare_truth_real_widebins_dcaz.C'("'${particle_array[iparticle]}'",'${etabin_array[i]}','${etabin_array[i+1]}')'
root -b -l -q doCompare_truth_real_widebins_dcaT.C'("'${particle_array[iparticle]}'",'${etabin_array[i]}','${etabin_array[i+1]}')'
done
done
