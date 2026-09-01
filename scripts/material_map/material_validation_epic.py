#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024 Shujie Li

import os
import argparse
from pathlib import Path

import acts
import acts.examples.dd4hep
from acts.examples import Sequencer

import epic
from material_validation import runMaterialValidation

u = acts.UnitConstants

if "__main__" == __name__:

    p = argparse.ArgumentParser(
        description="Script to produce propogation validation for ePIC material mapping."
    )
    p.add_argument(
        "--xmlFile",
        default=os.environ.get("DETECTOR_PATH", "") + os.environ.get("DETECTOR_CONFIG", "") + ".xml",
        help="input xml file containing ePIC geometry",
    )
    p.add_argument(
        "--matFile",
        type=str,
        default="",
        help="input material map file with extension, can be either xx.json or xx.cbor (optional)",
    )
    p.add_argument(
        "--outputName",
        type=str,
        default="propagation-material",
        help="customized name of the output rootfile (without .root extension)",
    )
    p.add_argument(
        "-n","--nevents",
        type=int,
        default=1000,
        help="number of events to run",
    )

    p.add_argument(
        "-t","--ntracks",
        type=int,
        default=1000,
        help="number of tracks per event")

    args = p.parse_args()

    detector = epic.getDetector(args.xmlFile, args.matFile)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    s = Sequencer(events=args.nevents, numThreads=-1)

    runMaterialValidation(
        surfaces=materialSurfaces,
        s=s,
        tracksPerEvent=args.ntracks,
        outputFileBase=os.path.join(os.getcwd(), args.outputName),
        trackingGeometry=trackingGeometry,
    ).run()
