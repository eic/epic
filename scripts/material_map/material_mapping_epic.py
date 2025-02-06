#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024 Shujie Li

import os
import argparse

import acts
from acts.examples import JsonFormat

import epic
from material_mapping import runMaterialMapping

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

    detector, trackingGeometry, decorators = epic.getDetector(
        args.xmlFile, args.geoFile)


    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir = os.getcwd(),
        inputDir  = os.getcwd(),
        readCachedSurfaceInformation=False,
        mapVolume= False,
        mapName  = mapName,
        mapFormat = mapFormat,
    ).run()
