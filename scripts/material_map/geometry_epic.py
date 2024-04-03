#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024 Shujie Li

import os
import argparse

import acts

import epic
from geometry import runGeometry

if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Script to generate geometry-map.json for ePIC geometry"
    )
    p.add_argument(
        "-i",
        "--xmlFile",
        default=(
            os.environ.get("DETECTOR_PATH", "")
            + "/"
            + os.environ.get("DETECTOR_CONFIG", "")
            + ".xml"
        ),
        help="Input xml file containing ePIC geometry",
    )
    args = p.parse_args()

    detector, trackingGeometry, decorators = epic.getDetector(args.xmlFile)

    runGeometry(
        trackingGeometry,
        decorators,
        outputDir=os.getcwd(),
        outputObj=False,
        outputCsv=False,
        outputJson=True,
        outputRoot=True,
    )
