#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024 Shujie Li

import os
import argparse
from pathlib import Path

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
    p.add_argument(
        "-o",
        "--outputName",
        type=str,
        default="geometry-map",
        help="output name for the geometry JSON file (without .json extension)",
    )
    args = p.parse_args()

    detector = epic.getDetector(args.xmlFile)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    import tempfile
    import shutil

    # Generate geometry to a temporary directory, then rename the output file
    tmpdir = Path(tempfile.mkdtemp())
    try:
        runGeometry(
            trackingGeometry,
            decorators,
            outputDir=tmpdir,
            outputObj=False,
            outputCsv=False,
            outputSurfacesJson=True,
        )
        # Move the generated geometry-map.json to the desired output name
        src_file = tmpdir / "geometry-map.json"
        dst_file = Path.cwd() / f"{args.outputName}.json"
        if src_file.exists():
            shutil.move(str(src_file), str(dst_file))
        else:
            raise FileNotFoundError(f"Expected geometry file not found at {src_file}")
    finally:
        shutil.rmtree(tmpdir)
