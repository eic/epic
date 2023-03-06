# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Chao Peng
'''
    A script to scan the materials and draw a plot for it
    It reads a json configuration file that defines the scan range, step, and the material-detector correspondence
'''

import os
import json
import argparse
import subprocess
import dd4hep
import DDRec
import pandas as pd
import numpy as np
from io import StringIO
from collections import OrderedDict


'''
    re-implementation of MaterialScan::Print in DD4Hep::DDRec
    MaterialScan does not have a python interface and that function is relatively simple, so here it is
    mng: instance of DDRec::MaterialManager
    start: 3D vector for start point
    end: 3D vector for end point
    epsilon: step size
'''
def material_scan(mng, start, end, epsilon=1e-4):
    p0 = np.array(start)
    p1 = np.array(end)
    direction = (p1 - p0)/np.linalg.norm(p1 - p0)
    placements = mng.placementsBetween(tuple(p0), tuple(p1), epsilon);

    # calculate material layer by layer
    int_x0 = 0
    int_lambda = 0
    path_length = 0
    res = []
    cols = [
        'material', 'Z', 'A', 'density',
        'rad_length', 'int_length', 'thickness', 'path_length',
        'int_X0', 'int_lamda', 'end_x', 'end_y', 'end_z'
        ]
    for pv, l in placements:
        # print(pv, l)
        path_length += l
        pcurr = p0 + path_length*direction
        mat = pv.GetMedium().GetMaterial()
        radl = mat.GetRadLen()
        intl = mat.GetIntLen()
        int_x0 += l/radl
        int_lambda += l/intl
        res.append([
            mat.GetName(), mat.GetZ(), mat.GetA(), mat.GetDensity(), radl, intl,
            l, path_length, int_x0, int_lambda, pcurr[0], pcurr[1], pcurr[2]
            ])
    dft = pd.DataFrame(data=res, columns=cols)
    return dft


'''
    A parser function to convert output from materialScan to a pandas dataframe
    compact: is path to compact file
    start_point: a list or tuple for 3D coordinates (cm, according to materialScan)
    end_point: a lit or tuple for 3D coordinates (cm, according to materialScan)
    return: a dataframe for materialScan results
'''
# NOTE: it is much slower than material_scan, so never used
def material_scan2(compact, start_point, end_point, timeout=200):
    EXEC_NAME = 'materialScan'
    cmd = '{} {} {} {} {} {} {} {}'.format(EXEC_NAME, compact, *np.hstack([start_point, end_point]))
    byteout = subprocess.check_output(cmd.split(' '), timeout=timeout)
    lines = []

    # NOTE: the code below is for materialScan output as of 03/05/2023
    # it may need change if the output format is changed
    first_digit = False
    for i, l in enumerate(byteout.decode('UTF-8').rstrip().split('\n')):
        line = l.strip('| ')
        if not first_digit and not line[:1].isdigit():
            continue
        first_digit = True
        if not line[:1].isdigit():
            break
        # break coordinates for endpoint, which has a format of (x, y, z)
        lines.append(line.strip('| ').translate({ord(i): None for i in '()'}).replace(',', ' '))

    cols = [
        'material', 'Z', 'A', 'density',
        'rad_length', 'int_length', 'thickness', 'path_length',
        'int_X0', 'int_lamda', 'end_x', 'end_y', 'end_z'
        ]

    dft = pd.read_csv(StringIO('\n'.join(lines)), sep='\s+', header=None, index_col=0, names=cols)
    return dft.astype({key: np.float64 for key in cols[1:]})


DEFAULT_CONFIG = dict(
    eta_range=[-1.5, 1.5, 0.5],
    phi=0.,
    path_r=118,
    start_point=[0., 0., 0.],
    detectors=OrderedDict([
        ('BECal ScFi', dict(
            materials=['SciFiPb*'],
            geocuts=[
                'math.sqrt(x**2 + y**2) >= {EcalBarrel_rmin}',
                'math.sqrt(x**2 + y**2) <= {EcalBarrel_SensitiveLayers_rmax}',
                'abs(z - {EcalBarrel_Calorimeter_offset}) <= {EcalBarrel_Calorimeter_length}/2.'],
            )),
        ])
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            dest='compact',
            help='Top-level xml file of the detector description.'
            )
    parser.add_argument(
            '-c', '--config', default=None,
            help='A JSON configuration file.'
            )
    args = parser.parse_args()

    if not os.path.exists(args.compact):
        print('Cannot find {}'.format(args.compact))
        exit(-1)

    # get configurations
    config = DEFAULT_CONFIG.copy()
    if args.config is not None:
        with open(args.config) as json_file:
            config = json.load(args.config)
    if 'detectors' not in config:
        print('No detectors defined in configuration file, draw all materials')
        config['detectors'] = OrderedDict([
            ('All Materials', dict(materials=['*'], geocuts=[])),
            ])

    # geometry initialization
    desc = dd4hep.Detector.getInstance()
    desc.fromXML(args.compact)
    mat_manager = DDRec.MaterialManager(desc.worldVolume())

    # build a constants dictionary for geocuts
    det_constants = {}
    for key, _ in desc.constants():
        try:
            det_constants[key] = desc.constantAsDouble(key)
        except Exception:
            det_constants[key] = desc.constantAsString(key)

    # update and compile geocuts
    # NOTE: variables below (x, y, z) will be used by compiled cuts
    x, y, z = 0., 0., 0.
    for det, dconf in config['detectors'].items():
        geocuts = []
        for gcut in dconf['geocuts']:
            cutstr = gcut.format(**det_constants)
            # print(cutstr)
            geocuts.append(compile(cutstr, '<string>', 'eval'))
        config['detectors'][det]['geocuts'] = geocuts

    # scan materials
    eta_range = config.get('eta_range')
    etas = np.arange(*eta_range)
    if etas[-1] < eta_range[1]:
        etas = np.hstack([etas, eta_range[1]])
    r = config.get('path_r')
    phi = config.get('phi')
    start_point = config.get('start_point')
    for eta in etas:
        x = r*np.cos(phi/180.*np.pi)
        y = r*np.sin(phi/180.*np.pi)
        z = r*np.sinh(eta)
        # print('({:.2f}, {:.2f}, {:.2f})'.format(x, y, z))
        print(material_scan(mat_manager, start_point, (x, y, z)).groupby('material')['thickness'].sum())

