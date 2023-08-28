#!/usr/bin/env python3

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Chao Peng
'''
    A script to use g4MaterialScan and collect the results
'''

import os
import argparse
import subprocess
import pandas as pd
import numpy as np
from io import StringIO
from collections import OrderedDict as odict

pd.set_option('display.max_rows', 1000)
PROGRESS_STEP = 10


'''
    A parser function to convert output from g4MaterialScan to a pandas dataframe
    compact: is path to compact file
    start_point: a list or tuple for 3D coordinates (e.g., 0,0,0)
    direction: a list or tuple for direction (e.g., 0.1,0.2,1.0)
    return: a dataframe for materialScan results
'''
def g4_material_scan(compact, start_point, direction, timeout=200):
    EXEC_NAME = 'g4MaterialScan'
    cmd = '{} --compact={} --position=\"{},{},{}\" --direction=\"{},{},{}\"'\
          .format(EXEC_NAME, compact, *np.hstack([start_point, direction]))
    output = os.popen(cmd).read()
    # output = subprocess.check_output(cmd.split(' '), timeout=timeout).decode('UTF-8')

    # find material scan lines
    lines = []
    add_line = False

    for l in output.split('\n'):
        if add_line:
            lines.append(l.strip())
        if l.strip().startswith('MaterialScan'):
            add_line = True


    # NOTE: the code below is for materialScan output as of 03/05/2023
    # it may need change if the output format is changed
    scans = []
    first_digit = False
    for i, l in enumerate(lines):
        line = l.strip('| ')
        if not first_digit and not line[:1].isdigit():
            continue
        first_digit = True
        if not line[:1].isdigit():
            break
        # break coordinates for endpoint, which has a format of (x, y, z)
        scans.append(line.strip('| ').translate({ord(i): None for i in '()'}).replace(',', ' '))

    cols = [
        'material', 'Z', 'A', 'density',
        'rad_length', 'int_length', 'thickness', 'path_length',
        'int_X0', 'int_lamda', 'end_x', 'end_y', 'end_z'
        ]

    dft = pd.read_csv(StringIO('\n'.join(scans)), sep='\s+', header=None, index_col=0, names=cols)
    return dft.astype({key: np.float64 for key in cols[1:]})


'''
    A helper function to convert a string (<min>[:<max>[:<step>]]) to an array
'''
def args_array(arg, step=1, include_end=True):
    vals = [float(x.strip()) for x in arg.split(':')]
    # empty or only one value
    if len(vals) < 2:
        return np.array(vals)
    # has step input
    if len(vals) > 2:
        step = vals[2]
    # inclusion of the endpoint (max)
    if include_end:
        vals[1] += step
    return np.arange(vals[0], vals[1], step)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            dest='compact',
            help='Top-level xml file of the detector description.'
            )
    parser.add_argument(
            '-o', '--output', default='g4_materials.csv',
            help='A csv file contains the scan info.'
            )
    parser.add_argument(
            '--start-point', default='0,0,0',
            help='Start point of the scan, use the format \"x,y,z\", unit is cm.'
            )
    parser.add_argument(
            '--eta', default='-4.0:4.0:0.1',
            help='Eta range, in the format of \"<min>[:<max>[:<step>]]\".'
            )
    parser.add_argument(
            '--phi', default='0:30:1',
            help='Phi angle range, in the format of \"<min>[:<max>[:<step>]]\" (degree).'
            )
    args = parser.parse_args()

    if not os.path.exists(args.compact):
        print('Cannot find {}'.format(args.compact))
        exit(-1)

    start_point = np.array([float(v.strip()) for v in args.start_point.split(',')])
    etas = args_array(args.eta)
    phis = args_array(args.phi)
    # sanity check
    if not len(phis):
        print('No phi values from the input {}, aborted!'.format(args.phi))
        exit(-1)
    mats_indices = odict()
    # a data buffer for the X0 values of (eta, material), 50 should be large enough
    data = np.zeros(shape=(len(etas), 50))
    # scan eta
    for i, eta in enumerate(etas):
        if i % PROGRESS_STEP == 0:
            print('Scanned {:d}/{:d} lines for {:.2f} < eta < {:.2f}, each line contains {:d} phi angles from {:.2f} to {:.2f}'\
                  .format(i, len(etas), etas[0], etas[-1], len(phis), phis[0], phis[-1])
                  )
        # average over phi
        eta_scan = pd.DataFrame()
        for phi in phis:
            direction = (np.cos(phi/180.*np.pi), np.sin(phi/180.*np.pi), np.sinh(eta))
            dfa = g4_material_scan(args.compact, start_point, direction)
            dfa.loc[:, 'X0'] = dfa['int_X0'].diff(1).fillna(dfa['int_X0'])
            phi_scan = dfa.groupby('material')['X0'].sum().to_frame(name=phi)
            # using pd.DataFrame.merge to combine results, using outer so new materials at specific phi angle can be added in
            eta_scan = eta_scan.merge(phi_scan, how='outer', left_index=True, right_index=True)
        # print(eta_scan)
        # print(eta_scan.sum(axis=1))
        # replace nan value with 0 and then calculate average values
        for mat, xval in eta_scan.fillna(0.).mean(axis=1).items():
            if mat not in mats_indices:
                mats_indices[mat] = len(mats_indices)
            j = mats_indices.get(mat)
            data[i, j] = xval

    result = pd.DataFrame(columns=mats_indices.keys(), index=etas, data=data[:, :len(mats_indices)])
    result.to_csv(args.output)
