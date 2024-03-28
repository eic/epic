#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024 Shujie Li

## Stand alone function to build ePIC geometry with ACTS python bindings
## for material mapping
## Shujie Li, 03, 2024

from pathlib import Path

import acts
import acts.examples.dd4hep

from acts import (
    Vector4,
    MaterialMapJsonConverter
)

import json

def getDetector(
    xmlFile,
    jsonFile="",
    logLevel=acts.logging.WARNING,
):
    customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)
    logger = acts.logging.getLogger("epic.getDetector")

    matDeco = None
    if len(jsonFile)>0:
        file = Path(jsonFile)
        logger.info("Adding material from %s", file.absolute())
        matDeco = acts.IMaterialDecorator.fromFile(
            file,
            level=customLogLevel(maxLevel=acts.logging.INFO),
        )

    dd4hepConfig = acts.examples.dd4hep.DD4hepGeometryService.Config(
        xmlFileNames=[xmlFile],
        logLevel=logLevel,
        dd4hepLogLevel=customLogLevel(),
    )
    detector = acts.examples.dd4hep.DD4hepDetector()

    config = acts.MaterialMapJsonConverter.Config()

    trackingGeometry, deco = detector.finalize(dd4hepConfig, matDeco)

    return detector, trackingGeometry, deco
