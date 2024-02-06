PTMOMLOW=$1
PTMOMHIGH=$2
THETALOW=$3
THETAHIGH=$4
DETXML=$5
SLCTPART=$6
#ddsim --runType vis --compactFile $DETECTOR_PATH/epic_${DETXML}.xml --macro macro/vis.mac --outputFile test.edm4hep.root -G --gun.thetaMin "$THETALOW*deg" --gun.thetaMax "$THETAHIGH*deg" --gun.distribution "uniform" --gun.momentumMin "$PTMOMLOW*GeV" --gun.momentumMax "$PTMOMHIGH*GeV" --gun.particle "$SLCTPART" -N 1

ddsim --runType vis --compactFile $DETECTOR_PATH/epic_${DETXML}.xml --macro macro/vis.mac --outputFile test.edm4hep.root -G --gun.thetaMin "$THETALOW*deg" --gun.thetaMax "$THETAHIGH*deg" --gun.distribution "pseudorapidity" --gun.momentumMin "$PTMOMLOW*GeV" --gun.momentumMax "$PTMOMHIGH*GeV" --gun.particle "$SLCTPART" -N 0

#ddsim --runType vis --compactFile $DETECTOR_PATH/epic_${DETXML}.xml --macro macro/vis.mac --outputFile test.edm4hep.root --numberOfEvents 1 --inputFiles /home/niviths/Downloads/pythia8NCDIS_10x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv.hepmc

#-numberOfEvents 25 --inputFiles input.hepmc  --gun.phiMin 1.58 --gun.phiMax 1.6

#e.g.
#bash plot_geo_wGun.sh 1.0 50.0 4. 25. gfhcal_only pi-
#bash plot_geo_wGun.sh 20.0 20.0 17. 18. gfhcal_only pi-

#bash plot_geo_wGun.sh 1.0 30.0 4. 25. arches mu-

#bash plot_geo.sh 0.1 20.0 4. 176. arches_trks_wTOF pi-
#--enableG4GPS

# sim new event
#/run/beamOn

#qt or vis

#bash plot_geo_wGun.sh 20.0 20.0 14. 15. lfhcal_testbeam pi-
#bash plot_geo_wGun.sh 1.0 50.0 4. 25. tof_endcap_only pi-
#bash plot_geo_wGun.sh 1.0 50.0 4. 25. craterlake pi-
