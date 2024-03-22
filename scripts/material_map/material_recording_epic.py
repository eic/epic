#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024 Shujie Li

import os
import warnings
from pathlib import Path
import argparse

import acts
from acts.examples import (
    GaussianVertexGenerator,
    ParametricParticleGenerator,
    FixedMultiplicityGenerator,
    EventGenerator,
    RandomNumbers,
)

import acts.examples.dd4hep
import acts.examples.geant4
import acts.examples.geant4.dd4hep

import epic
from material_recording import runMaterialRecording

u = acts.UnitConstants

_material_recording_executed = False


def main():

    p = argparse.ArgumentParser()
    p.add_argument(
        "-n", "--events", type=int, default=1000, help="Number of events to generate"
    )
    p.add_argument(
        "-t", "--tracks", type=int, default=100, help="Particle tracks per event"
    )
    p.add_argument(
        "-i", "--xmlFile", type=str, default=os.environ.get("DETECTOR_PATH", "") + os.environ.get("DETECTOR_CONFIG", "") + ".xml", help="DD4hep input file"
    )
    p.add_argument(
        "-o", "--outputName", type=str, default="geant4_material_tracks.root", help="Name of the output rootfile"
    )
    p.add_argument(
        "--eta_min",
        type=float,
        default=-8.0,
        help="eta min (optional)",
    )
    p.add_argument(
        "--eta_max",
        type=float,
        default=8.0,
        help="eta max (optional)",
    )
    args = p.parse_args()

    detector, trackingGeometry, decorators = epic.getDetector(
        args.xmlFile)

    detectorConstructionFactory = (
        acts.examples.geant4.dd4hep.DDG4DetectorConstructionFactory(detector)
    )

    runMaterialRecording(
        detectorConstructionFactory=detectorConstructionFactory,
        tracksPerEvent=args.tracks,
        outputDir=os.getcwd(),
        etaRange=(args.eta_min, args.eta_max),
        s=acts.examples.Sequencer(events=args.events, numThreads=1),
    ).run()


if "__main__" == __name__:
    main()
