#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024 Shujie Li

import os
import argparse
from pathlib import Path

import acts
from material_mapping import runMaterialMapping

try:
    from acts.examples.json import JsonFormat
except ImportError:
    from acts.examples import JsonFormat

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
        "--matFileBase",
        type=str,
        default="material-map",
        help="base name for the generated material map (without extension)",
    )
    p.add_argument(
        "--matFileFormat",
        type=str,
        nargs="+",
        default=["json", "root"],
        help="output format(s) for the material map (json, cbor, root). Default: json, root",
    )
    p.add_argument(
        "--inputRootFile",
        type=str,
        default="geant4_material_tracks.root",
        help="input ROOT file with material tracks (typically from material_recording_epic.py)",
    )
    args = p.parse_args()

    # Parse material map formats
    mapFormats = args.matFileFormat
    # Validate formats
    valid_formats = {"json", "cbor", "root"}
    for fmt in mapFormats:
        if fmt not in valid_formats:
            print(f'ERROR(material_mapping_epic.py): invalid format "{fmt}". Must be one of: {", ".join(valid_formats)}')
            exit(1)

    detector = epic.getDetector(args.xmlFile)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    outputFileBase = os.path.join(os.getcwd(), args.matFileBase)

    runMaterialMapping(
        surfaces=materialSurfaces,
        inputFile=Path(args.inputRootFile),
        outputFileBase=outputFileBase,
        outputMapFormats=mapFormats,
    ).run()
