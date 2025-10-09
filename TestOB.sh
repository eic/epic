#!/bin/bash
ddsim --compactFile epic_craterlake_tracking_only_cut.xml --outputFile test.edm4hep.root --numberOfEvents 1000 --enableGun --gun.thetaMin 50*deg --gun.thetaMax 130*deg --gun.distribution eta --gun.particle pi- --gun.momentumMin 10*GeV --gun.momentumMax 10*GeV --gun.multiplicity 1
eicrecon -Pnthreads=1 -Pjana:debug_plugin_loading=1 -Ppodio:output_file=reconTest.edm4eic.root -Pdd4hep:xml_files=epic_craterlake_tracking_only_cut.xml -Ppodio:output_collections="MCParticles,CentralCKFTrajectories,CentralCKFTrackParameters,CentralCKFSeededTrackParameters,CentralCKFTruthSeededTrackParameters,CentralTrackVertices,ReconstructedChargedParticles,ReconstructedChargedParticleAssociations" test.edm4hep.root
root -b TestOB.C reconTest.edm4eic.root
