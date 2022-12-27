#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Wouter Deconinck

from __future__ import absolute_import, unicode_literals
import os
import time
import logging

import argparse
parser = argparse.ArgumentParser(
     prog='checkOverlaps.py',
     description='''Check for overlaps using Geant4''',
     epilog='''
     This program checks the compact detector file for overlaps using Geant4.
         ''')
parser.add_argument("-c", "--compact", help="compact detector file",default="athena.xml")
parser.add_argument("-r", "--resolution", help="number of points on surface",default="10000")
parser.add_argument("-t", "--tolerance", help="minimum distance (in mm) to report overlaps",default="0.1")
parser.add_argument("-v", "--verbose", help="print output", action='store_true')

args = parser.parse_args()

import DDG4
from g4units import keV, GeV, mm, ns, MeV

def run():
  kernel = DDG4.Kernel()
  description = kernel.detectorDescription()
  kernel.loadGeometry(str("file:" + args.compact))

  DDG4.importConstants(description)

  geant4 = DDG4.Geant4(kernel)
  ui = geant4.setupCshUI(ui=None)
  ui.Commands = [
      '/geometry/test/resolution {}'.format(args.resolution),
      '/geometry/test/tolerance {}'.format(args.tolerance),
      '/geometry/test/verbosity {}'.format(1 if args.verbose else 0),
      '/geometry/test/run'
      ]
  kernel.configure()
  kernel.initialize()
  kernel.run()
  kernel.terminate()


if __name__ == "__main__":
  run()
