#!/usr/bin/env python
"""
DD4hep simulation with some argument parsing
Based on M. Frank and F. Gaede runSim.py
   @author  A.Sailer
   @version 0.1
Modified with standard EIC EPIC requirements.
"""
from __future__ import absolute_import, unicode_literals
import logging
import sys
import math

#from DDSim import SystemOfUnits
from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV, radian

SIM = DD4hepSimulation()
################################################################
# The compact XML file
# This defines the geometry version to use -->
# If making changes, you NEED to compile a local version
# of the ePIC repo!
#
#
# FIXME: You MUST source the right setup file is you compile 
#        ePIC locally!
#
#        source epic/install/setup.sh
#
#
#        The default magnetic field settings are for 18x275 GeV
#################################################################

#SIM.compactFile = "/epic_ip6.xml"
## Lorentz boost for the crossing angle, in radian!

#############################################################
#   Particle Gun simulations use these options
#   --> comment out if using MC input
#############################################################

#SIM.enableDetailedShowerMode = False
SIM.enableG4GPS = False
SIM.enableG4Gun = False
SIM.enableGun = True
#SIM.gun.direction = (math.sin(-0.025*radian), 0, math.cos(-0.025*radian))
SIM.crossingAngleBoost = -0.025*radian


## InputFiles for simulation .stdhep, .slcio, .HEPEvt, .hepevt, .hepmc, .pairs files are supported
SIM.inputFiles = []
## Macro file to execute for runType 'run' or 'vis'
SIM.macroFile = ""
## number of events to simulate, used in batch mode
SIM.numberOfEvents = 5000
## Outputfile from the simulation,only lcio output is supported
SIM.outputFile = "RP_protons_testing_275_GeV_5k_events_8_29_2023_0.edm4hep.root"
## Physics list to use in simulation
#SIM.physicsList = None
## Verbosity use integers from 1(most) to 7(least) verbose
## or strings: VERBOSE, DEBUG, INFO, WARNING, ERROR, FATAL, ALWAYS
#SIM.printLevel = "VERBOSE"
## The type of action to do in this invocation
## batch: just simulate some events, needs numberOfEvents, and input file or gun
## vis: enable visualisation, run the macroFile if it is set
## run: run the macroFile and exit
## shell: enable interactive session
SIM.runType = "run"
## Skip first N events when reading a file
SIM.skipNEvents = 0
## Steering file to change default behaviour
#SIM.steeringFile = None
## FourVector of translation for the Smearing of the Vertex position: x y z t
SIM.vertexOffset = [0.0, 0.0, 0.0, 0.0]
## FourVector of the Sigma for the Smearing of the Vertex position: x y z t
SIM.vertexSigma = [0.0, 0.0, 0.0, 0.0]


################################################################################
## Action holding sensitive detector actions
##   The default tracker and calorimeter actions can be set with
##
##   >>> SIM = DD4hepSimulation()
##   >>> SIM.action.tracker = ('Geant4TrackerWeightedAction', {'HitPositionCombination': 2, 'CollectSingleDeposits': False})
##   >>> SIM.action.calo    = "Geant4CalorimeterAction"
##
##   for specific subdetectors specific sensitive detectors can be set based on pattern matching
##
##   >>> SIM = DD4hepSimulation()
##   >>> SIM.action.mapActions['tpc'] = "TPCSDAction"
##
##   and additional parameters for the sensitive detectors can be set when the map is given a tuple
##
##   >>> SIM = DD4hepSimulation()
##   >>> SIM.action.mapActions['ecal'] =( "CaloPreShowerSDAction", {"FirstLayerNumber": 1} )
##
##
################################################################################

##  set the default calorimeter action
##SIM.action.calo = "Geant4ScintillatorCalorimeterAction"

##  create a map of patterns and actions to be applied to sensitive detectors
##         example: SIM.action.mapActions['tpc'] = "TPCSDAction"
##SIM.action.mapActions = {}

##SIM.action.mapActions['tpc'] = "TPCSDAction"

##  set the default tracker action
##SIM.action.tracker = ('Geant4TrackerWeightedAction', {'HitPositionCombination': 2, 'CollectSingleDeposits': False})


################################################################################
## Configuration for the magnetic field (stepper)
################################################################################
SIM.field.delta_chord = 1e-03
SIM.field.delta_intersection = 1e-03
SIM.field.delta_one_step = .5e-01*mm
SIM.field.eps_max = 1e-03
SIM.field.eps_min = 1e-04
#SIM.field.equation = "Mag_UsualEqRhs"
SIM.field.largest_step = 100.0*mm
SIM.field.min_chord_step = 1.e-2*mm
SIM.field.stepper = "HelixSimpleRunge"


################################################################################
## Configuration for sensitive detector filters
##
##   Set the default filter for tracker or caliromter
##   >>> SIM.filter.tracker = "edep1kev"
##   >>> SIM.filter.calo = ""
##
##   Assign a filter to a sensitive detector via pattern matching
##   >>> SIM.filter.mapDetFilter['FTD'] = "edep1kev"
##
##   Or more than one filter:
##   >>> SIM.filter.mapDetFilter['FTD'] = ["edep1kev", "geantino"]
##
##   Don't use the default filter or anything else:
##   >>> SIM.filter.mapDetFilter['TPC'] = None ## or "" or []
##
##   Create a custom filter. The dictionary is used to instantiate the filter later on
##   >>> SIM.filter.filters['edep3kev'] = dict(name="EnergyDepositMinimumCut/3keV", parameter={"Cut": 3.0*keV} )
##
##
################################################################################

##  default filter for calorimeter sensitive detectors; this is applied if no other filter is used for a calorimeter
##SIM.filter.calo = "edep0"

##  list of filter objects: map between name and parameter dictionary
#SIM.filter.filters[ 'edep100MeV' ] = dict( name="EnergyDepositMinimumCut/100MeV", parameter={"Cut" : 1000*MeV} )
#SIM.filter.filters = {'alexEDep': {'parameter': {'Cut': 100.0}, 'name': 'EnergyDepositMinimumCut/cut0'}}

#SIM.filter.calo = "edep100MeV"
#SIM.filter.mapDetFilter['ForwardRomanPot_Station_1'] = "edep100MeV"
#SIM.filter.mapDetFilter['ForwardRomanPot_Station_2'] = "edep100MeV"
#SIM.part.minimalKineticEnergy = 100*MeV 

##  a map between patterns and filter objects, using patterns to attach filters to sensitive detector
##SIM.filter.mapDetFilter = {}

##SIM.filter.mapDetFilter['TPC'] = "edep0"

##  default filter for tracking sensitive detectors; this is applied if no other filter is used for a tracker
##SIM.filter.tracker = "edep1kev"


################################################################################
## Configuration for the GuineaPig InputFiles
################################################################################

## Set the number of pair particles to simulate per event.
##     Only used if inputFile ends with ".pairs"
##     If "-1" all particles will be simulated in a single event
##
#SIM.guineapig.particlesPerEvent = "1"


################################################################################
## Configuration for the DDG4 ParticleGun
################################################################################

##  direction of the particle gun, 3 vector
##  shoot in the barrel region
#SIM.gun.direction = (0, 1, 0)

## choose the distribution of the random direction for theta
##
##     Options for random distributions:
##
##     'uniform' is the default distribution, flat in theta
##     'cos(theta)' is flat in cos(theta)
##     'eta', or 'pseudorapidity' is flat in pseudorapity
##     'ffbar' is distributed according to 1+cos^2(theta)
##
##     Setting a distribution will set isotrop = True
##
SIM.gun.distribution = 'uniform'
SIM.gun.energy = 275.0*GeV ## default energy value is MeV
#SIM.gun.momentumMin = 80*GeV
#SIM.gun.momentumMax = 100*GeV

##  isotropic distribution for the particle gun
##
##     use the options phiMin, phiMax, thetaMin, and thetaMax to limit the range of randomly distributed directions
##     if one of these options is not None the random distribution will be set to True and cannot be turned off!
##
SIM.gun.isotrop = True
#SIM.gun.multiplicity = 1
SIM.gun.particle = "proton"
##SIM.gun.phiMax = pi

## Minimal azimuthal angle for random distribution
##SIM.gun.phiMin = 0

##  position of the particle gun, 3 vector. unit mm
##SIM.gun.position = (0, 0, 0)
SIM.gun.thetaMin = 0.002*radian
SIM.gun.thetaMax = 0.004*radian


################################################################################
## Configuration for the output levels of DDG4 components
################################################################################

## Output level for input sources
SIM.output.inputStage = 3

## Output level for Geant4 kernel
SIM.output.kernel = 3

## Output level for ParticleHandler
SIM.output.part = 3

## Output level for Random Number Generator setup
SIM.output.random = 6


################################################################################
## Configuration for the Particle Handler/ MCTruth treatment
################################################################################

## Enable lots of printout on simulated hits and MC-truth information
SIM.part.enableDetailedHitsAndParticleInfo = False

##  Keep all created particles
SIM.part.keepAllParticles = True

## Minimal distance between particle vertex and endpoint of parent after
##     which the vertexIsNotEndpointOfParent flag is set
##
SIM.part.minDistToParentVertex = 2.2e-14

## MinimalKineticEnergy to store particles created in the tracking region
#SIM.part.minimalKineticEnergy = 1*MeV

##  Printout at End of Tracking
SIM.part.printEndTracking = False

##  Printout at Start of Tracking
SIM.part.printStartTracking = False

## List of processes to save, on command line give as whitespace separated string in quotation marks
SIM.part.saveProcesses = ['Decay']


################################################################################
## Configuration for the PhysicsList
################################################################################
SIM.physics.decays = True
SIM.physics.list = "FTFP_BERT" # "FTFP_BERT"

##  location of particle.tbl file containing extra particles and their lifetime information
##


##  The global geant4 rangecut for secondary production
##
##     Default is 0.7 mm as is the case in geant4 10
##
##     To disable this plugin and be absolutely sure to use the Geant4 default range cut use "None"
##
##     Set printlevel to DEBUG to see a printout of all range cuts,
##     but this only works if range cut is not "None"
##
#SIM.physics.rangecut =  10.0*mm


################################################################################
## Properties for the random number generator
################################################################################

## If True, calculate random seed for each event based on eventID and runID
## allows reproducibility even when SkippingEvents
SIM.random.enableEventSeed = False
SIM.random.file = None
SIM.random.luxury = 1
SIM.random.replace_gRandom = True
SIM.random.seed = None
SIM.random.type = None 
