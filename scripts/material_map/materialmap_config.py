#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024 Shujie Li

## Generate ePIC material map for ACTS
## read config-map.json and turn on mapping for approach 1 and 2 of each sensitive surface.
import pandas as pd
import numpy as np
import json
import os
import argparse

if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Script to turn on all approach 1 and 2, and also the beampipe surface in config json file for matieral mapping"
    )
    p.add_argument(
        "-i","--inputFile",
        type=str,
        default="config-map.json",
        help=" input json file to be modified",
    )
    p.add_argument(
        "-o","--outputFile",
        type=str,
        default="config-map_new.json",
        help=" output json file",
    )

args     = p.parse_args()
print(args)
fname    = args.inputFile
out_name = args.outputFile


## load json file
f  = open(fname)
dd = json.load(f)

ee=dd['Volumes']['entries']

## print volume name and ID
print ("Volume ID        Name       Approaches")
for vv in np.arange(len(ee)):
    nn = ee[vv]['value']['NAME']

    if  "|" not in nn and "Gap"  not in nn:
        print(ee[vv]['volume'], nn,"1, 2")#print(ee[vv]['value'])#['NAME'])
        if "acts_beampipe_central::Barrel" in nn:
            v_beampipe = vv+1
            print(v_beampipe, nn, "X")

## find apporach 1 and 2 to turn on mapping
for vv in np.arange(1,1+len(dd['Surfaces'])):
    for ii,tt in enumerate(dd['Surfaces'][str(vv)]):
        if 'approach' in tt:
            dd['Surfaces'][str(vv)][ii]['value']['material']['mapMaterial']=True
        ## turn on beampipe surface and defind binning
        elif vv==v_beampipe:
            if tt['value']['bounds']['type']=='CylinderBounds':
            # print (dd['Surfaces'][str(vv)][ii])
                dd['Surfaces'][str(vv)][ii]['value']['material']['mapMaterial']=True
                dd['Surfaces'][str(vv)][ii]['value']['material']['binUtility']['binningdata'][0]['bins']=36
                dd['Surfaces'][str(vv)][ii]['value']['material']['binUtility']['binningdata'][1]['bins']=200


with open(out_name, "w") as outfile:
    json.dump(dd, outfile, indent=4)

print("Done! Updated config file at "+out_name)
