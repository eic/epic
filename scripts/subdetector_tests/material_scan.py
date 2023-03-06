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
import pandas as pd
import numpy as np
from io import StringIO
from collections import OrderedDict
from matplotlib import pyplot as plt

import ROOT
import dd4hep
import DDRec

'''
    re-implementation of MaterialScan::Print in DD4Hep::DDRec
    MaterialScan does not have a python interface and that function is relatively simple, so here it is
    desc: DD4hep::Detector
    start: 3D vector for start point
    end: 3D vector for end point
    epsilon: step size
'''
def material_scan(desc, start, end, epsilon=1e-5):
    rvec = ROOT.Math.XYZVector()
    mat_mng = DDRec.MaterialManager(desc.worldVolume())
    # id_conv = DDRec.CellIDPositionConverter(desc)
    # det_dict = {d.id(): n for n, d in desc.world().children()}
    # print(det_dict)

    p0 = np.array(start)
    p1 = np.array(end)
    direction = (p1 - p0)/np.linalg.norm(p1 - p0)
    placements = mat_mng.placementsBetween(tuple(p0), tuple(p1), epsilon);

    # calculate material layer by layer
    int_x0 = 0
    int_lambda = 0
    path_length = 0
    res = []
    cols = [
        # 'detector',
        'material', 'Z', 'A', 'density',
        'radl', 'intl', 'thickness', 'path_length',
        'X0', 'lamda', 'end_x', 'end_y', 'end_z'
        ]
    for pv, l in placements:
        # print(pv, l)
        path_length += l
        pcurr = p0 + path_length*direction
        mat = pv.GetMedium().GetMaterial()
        radl = mat.GetRadLen()
        intl = mat.GetIntLen()
        x0 = l/radl
        lmd = l/intl
        # try to locate the detector
        # det_id = -1
        # try:
        #     rvec.SetCoordinates(pcurr)
        #     det_id = id_conv.findDetElement(rvec).id()
        # except Exception:
        #     pass
        res.append([
            # det_dict.get(det_id, 'Unknown'),
            mat.GetName(), mat.GetZ(), mat.GetA(), mat.GetDensity(), radl, intl,
            l, path_length, x0, lmd, pcurr[0], pcurr[1], pcurr[2],
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
# NOTE: it is much slower than material_scan, so only used for testing
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
    eta_range=[-2.0, 1.7, 0.01],
    phi=70.,
    path_r=120.,
    start_point=[0., 0., 0.],
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

    # geometry initialization
    desc = dd4hep.Detector.getInstance()
    desc.fromXML(args.compact)

    # scan materials
    eta_range = config.get('eta_range')
    etas = np.arange(*eta_range)
    path_r = config.get('path_r')
    phi = config.get('phi')
    start_point = config.get('start_point')
    dets = [
        # 'EcalBarrelImaging',
        'EcalBarrelScFi',
        ]
    detmats = ['SciFiPb_PbGlue', 'SciFiPb_Scintillator']
    # others for all detectors that are not listed above
    dets.append('Others')
    vals = np.zeros(shape=(len(etas), len(dets)))
    for i, eta in enumerate(etas):
        x = path_r*np.cos(phi/180.*np.pi)
        y = path_r*np.sin(phi/180.*np.pi)
        z = path_r*np.sinh(eta)
        # print('({:.2f}, {:.2f}, {:.2f})'.format(x, y, z))
        dfr = material_scan(desc, start_point, (x, y, z))
        dmap = {k: 'Others' for k in dfr['material'].unique() if k not in detmats}
        dmap.update({k: 'EcalBarrelScFi' for k in detmats})
        dfr.loc[:, 'det2'] = dfr['material'].map(dmap)
        x0_vals = dfr.groupby('det2')['X0'].sum().to_dict()
        for k, det in enumerate(dets):
            vals[i, k] = x0_vals.get(det, 0.)

    # plot
    fig, ax = plt.subplots(figsize=(16, 4), dpi=160, gridspec_kw={'top': 0.995, 'left': 0.08, 'right': 0.98})
    hbins = np.hstack([etas, etas[-1] + eta_range[2]]) - eta_range[2]/2.
    ax.hist(np.tile(etas, (len(dets), 1)).T, hbins, weights=vals, histtype='step',
            stacked=True, fill=True, alpha=1.0, label=dets, color=['royalblue', 'indianred'])
    ax.legend(loc="upper center", fontsize=22)
    ax.tick_params(direction='in', labelsize=22)
    ax.set_xlabel('$\eta$', fontsize=22)
    ax.set_ylabel('X0', fontsize=22)
    ax.grid(ls=':')
    ax.set_axisbelow(True)
    ax.set_xlim(-2.2, 1.8)
    # ax.set_yscale('log')
    fig.savefig('material_scan.png')
