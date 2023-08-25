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
    direction: a lit or tuple for direction (e.g., 0.1,0.2,1.0)
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
            '--eta-min', type=float, default=-4.0,
            help='Minimum eta for the scan.'
            )
    parser.add_argument(
            '--eta-max', type=float, default=4.0,
            help='Minimum eta for the scan.'
            )
    parser.add_argument(
            '--eta-nbins', type=int, default=801,
            help='Number of bins for the eta scan.'
            )
    parser.add_argument(
            '--phi', type=float, default=20.,
            help='Phi angle of the scan, unit is degree.'
            )
    args = parser.parse_args()

    if not os.path.exists(args.compact):
        print('Cannot find {}'.format(args.compact))
        exit(-1)

    start_point = np.array([float(v.strip()) for v in args.start_point.split(',')])
    etas = np.linspace(args.eta_min, args.eta_max, args.eta_nbins)
    mats_indices = odict()
    # a data buffer for the X0 values of (eta, material), 50 should be large enough
    data = np.zeros(shape=(len(etas), 50))
    for i, eta in enumerate(etas):
        if i % PROGRESS_STEP == 0:
            print('Scanned {:d}/{:d} lines for {:.2f} < eta < {:.2f}'.format(i, len(etas), etas[0], etas[-1]))
        direction = (np.cos(args.phi/180.*np.pi), np.sin(args.phi/180.*np.pi), np.sinh(eta))
        dfa = g4_material_scan(args.compact, start_point, direction)
        dfa.loc[:, 'X0'] = dfa['int_X0'].diff(1).fillna(dfa['int_X0'])
        for mat, xval in dfa.groupby('material')['X0'].sum().items():
            if mat not in mats_indices:
                mats_indices[mat] = len(mats_indices)
            j = mats_indices.get(mat)
            data[i, j] = xval

    result = pd.DataFrame(columns=mats_indices.keys(), index=etas, data=data[:, :len(mats_indices)])
    result.to_csv(args.output)
