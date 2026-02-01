#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024 Shujie Li

import os
import argparse
from pathlib import Path

import acts
from material_mapping import runMaterialMapping

try:
    from acts.examples.json import  JsonFormat
except ImportError:
    from acts.examples import  JsonFormat

import epic


if "__main__" == __name__:

    p = argparse.ArgumentParser(
        description="Script to generate material map for ePIC geometry"
    )
    p.add_argument(
        "--xmlFile",
        default=os.environ.get("DETECTOR_PATH", "")+"epic_craterlake.xml",
        help="input xml file containing ePIC geometry",
    )
    p.add_argument(
        "--geoFile",
        type=str,
        default="geometry-map.json",
        help="input json file to define volumes and layers used in material mapping",
    )
    p.add_argument(
        "--matFile",
        type=str,
        default="material-map.json",
        help="output filename for the generated material map, can be json and cbor formats",
    )
    args = p.parse_args()

    mapName = args.matFile.split('.')[0]
    if '.json' in args.matFile:
        mapFormat = JsonFormat.Json
    elif '.cbor' in args.matFile:
        mapFormat = JsonFormat.Cbor
    else:
        print('ERROR(material_mapping_epic.py): please provide a material map file in .json or .cbor format')
        exit()

    detector = epic.getDetector(
        args.xmlFile, args.geoFile)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir=Path.cwd(),
        inputDir=Path.cwd(),
        readCachedSurfaceInformation=False,
        mapVolume=False,
        mapName=mapName,
        mapFormat=mapFormat,
    ).run()
