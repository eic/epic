#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import absolute_import, unicode_literals
import os
import time
import logging

import argparse
parser = argparse.ArgumentParser(
     prog='convert_to_gdml.py',
     description='''Convert DD4Hep description to GDML''',
     epilog='''
     This program converts the compact detector file to a single GDML file.
         ''')
parser.add_argument("-c", "--compact", help="compact detector file",default="athena.xml")
parser.add_argument("-o", "--output", help="gdml detector file",default="athena.gdml")

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
  #
  # Setup the GDML writer action
  writer = DDG4.Action(kernel, 'Geant4GDMLWriteAction/Writer')
  writer.enableUI()
  kernel.registerGlobalAction(writer)
  ui.Commands = [
      '/ddg4/Writer/Output {}'.format(args.output),
      '/ddg4/Writer/OverWrite 1',
      '/ddg4/Writer/ModuleDepth 1',
      '/ddg4/Writer/write'
      ]
  kernel.configure()
  kernel.initialize()
  kernel.run()
  kernel.terminate()


if __name__ == "__main__":
  run()
