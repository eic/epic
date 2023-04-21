# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Chao Peng
'''
    A script to scan the materials and draw a plot for it
    It reads a json configuration file that defines the scan range, step, and the material-detector correspondence
'''

import os
import math
import json
import fnmatch
import argparse
import subprocess
import pandas as pd
import numpy as np
from io import StringIO
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.backends.backend_pdf

import ROOT
import dd4hep
import DDRec

pd.set_option('display.max_rows', 1000)
PROGRESS_STEP = 10


'''
    Re-implementation of MaterialScan::Print in DD4Hep::DDRec
    MaterialScan does not have a python interface and that function is relatively simple, so here it is
    desc: DD4hep::Detector
    start: 3D vector for start point
    end: 3D vector for end point
    epsilon: step size
'''
# NOTE: I tried to get detector information from placedVolume (or CellIDPositionConverter) but it is not useable
# (too many mistakes and wrong assignments).
def material_scan(desc, start, end, epsilon=1e-5):
    mat_mng = DDRec.MaterialManager(desc.worldVolume())
    # only use the top-level detectors
    dets = [d for n, d in desc.world().children()]
    # rvec = ROOT.Math.XYZVector()
    # id_conv = DDRec.CellIDPositionConverter(desc)
    # det_dict = {d.id(): n for n, d in desc.world().children()}

    p0 = np.array(start)
    p1 = np.array(end)
    direction = (p1 - p0)/np.linalg.norm(p1 - p0)
    print(p0, p1, direction)
    placements = mat_mng.placementsBetween(tuple(p0), tuple(p1), epsilon);

    # calculate material layer by layer
    int_x0 = 0
    int_lambda = 0
    path_length = 0
    # FIXME: this is a work-around to a known issue of dd4hep::rec::MaterialManager::placementsBetween
    # It is GEOMETRY DEPENDENT (the thickness of the first vacuum material layer)
    # It should be removed once the dd4hep issue is fixed, preparing a PR for that as of 04/21/2023
    # to check the thickness, see the difference between
    # materialScan epic_brycecanyon.xml 0 0 0 100  40  -0.01 | grep Vacuum
    # materialScan epic_brycecanyon.xml 0 0 0 100  40  0.01 | grep Vacuum
    if direction[2] < 0.:
        path_length += np.sqrt(1./(direction[0]**2 + direction[1]**2))*2.8
    # FIXME: work-around ended
    res = []
    for pv, l in placements:
        # print(pv, l)
        path_length += l
        pcurr = p0 + path_length*direction
        mat = pv.GetMedium().GetMaterial()
        radl = mat.GetRadLen()
        intl = mat.GetIntLen()
        x0 = l/radl
        lmd = l/intl
        d0 = 'Unknown'
        for d in dets:
            local = d.nominal().worldToLocal(pcurr)
            local = np.array([local.X(), local.Y(), local.Z()])
            if d.volume().Contains(local):
                d0 = d.GetName()
        # try to locate the detector
        # det_id = -1
        # try:
        #     rvec.SetCoordinates(pcurr)
        #     det_id = id_conv.findDetElement(rvec).id()
        # except Exception:
        #     pass
        res.append([
            d0,
            # det_dict.get(det_id, 'Unknown'),
            mat.GetName(), mat.GetZ(), mat.GetA(), mat.GetDensity(), radl, intl,
            l, path_length, x0, lmd,
            pcurr[0], pcurr[1], pcurr[2],
            # local[0], local[1], local[2],
            ])
    cols = [
        'detector',
        'material', 'Z', 'A', 'density',
        'radl', 'intl', 'thickness', 'path_length',
        'X0', 'lamda',
        'x', 'y', 'z', 'r_xy',
        # 'local_x', 'local_y', 'local_z'
        ]
    dft = pd.DataFrame(data=res, columns=cols)
    # print(dft[['detector', 'material', 'x', 'y', 'z', 'path_length', 'r_xy', 'local_x', 'local_y', 'local_z']].head(100))
    # print(dft.groupby('detector')['X0'].sum())
    return dft


# a default configuraiton for bareel ecal of epic_brycecanyon
DEFAULT_CONFIG = dict(
    eta_range=[-2.0, 1.7, 0.01],
    phi=20.,
    path_r=120.,
    start_point=[0., 0., 0.],
    # The order defines the priority for assigning material layers to a detector
    detectors=OrderedDict([
        ('EcalBarrelScFi', dict(
            # because of special materials, we don't need geo constraints
            materials=['SciFiPb*'],
            )),
        ('EcalBarrelImg', dict(
            materials=['Air', 'CarbonFiber', 'Silicon', 'Copper', 'Kapton', 'Epoxy'],
            geo=[
                'math.sqrt(x**2 + y**2) >= {EcalBarrel_rmin}',
                'math.sqrt(x**2 + y**2) <= {EcalBarrel_SensitiveLayers_rmax}',
                'abs(z - {EcalBarrel_Calorimeter_offset}) <= {EcalBarrel_Calorimeter_length}/2.',
                ],
            )),
        ('InnerTrackerSupport', dict(
            materials=['Aluminum', 'CarbonFiber'],
            geo=[
                'math.sqrt(x**2 + y**2) < {DIRC_rmin}',
                'z <= {TrackerSupportConeEndcapN_zmax}',
                'z >= -{TrackerSupportConeEndcapN_zmax}',
                ],
            )),
        # ('BarrelDIRC', dict(
        #     materials=['Aluminum', 'Quartz', 'Nlak33a'],
        #     geo=[
        #         'math.sqrt(x**2 + y**2) >= {DIRC_rmin}',
        #         'math.sqrt(x**2 + y**2) < {EcalBarrel_rmin}',
        #         'z <= {DIRCForward_zmax}',
        #         'z >= -{DIRCBackward_zmax}',
        #         ],
        #     )),
        # ('EcalEndcapN', dict(
        #     materials=['StainlessSteel', 'leadtungsten_optical', 'CarbonFiber', 'VM2000'],
        #     geo=[
        #         'math.sqrt(x**2 + y**2) < {EcalBarrel_rmin}',
        #         'z < -{EcalEndcapN_zmin}',
        #         ],
        #     )),
        # ('HcalEndcapN', dict(
        #     materials=['Polystyrene', 'Steel235'],
        #     geo=['z < {HcalEndcapN_zmin}'],
        #    )),
        # all other materials fall into this category
        ('Others', dict(
            materials=['*'],
            )),
        # end of detectors
        ])
    )

# execute the script
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

    # build a constants dictionary for geo constraints
    det_constants = {}
    for key, _ in desc.constants():
        try:
            det_constants[key] = desc.constantAsDouble(key)
        except Exception:
            det_constants[key] = desc.constantAsString(key)
    # print(json.dumps(det_constants, indent=4, sort_keys=True))

    # update and compile geo constraints
    # NOTE: variables below (x, y, z) will be used by compiled cuts
    x, y, z = 0., 0., 0.
    for det, dconf in config['detectors'].items():
        geo = []
        for g in dconf.get('geo', []):
            cutstr = g.format(**det_constants)
            # print('{} => {}'.format(g, cutstr))
            geo.append(compile(cutstr, '<string>', 'eval'))
        config['detectors'][det]['geo'] = geo

    # scan materials
    eta_range = config.get('eta_range')
    etas = np.arange(*eta_range)
    path_r = config.get('path_r')
    phi = config.get('phi')
    start_point = config.get('start_point')
    # dets = list(config.get('detectors').keys())
    dets = [n for n, _ in desc.world().children()]
    dets.append('Unknown')

    # array for detailed data
    # materials in configuration maybe using wild-card matching, so just assign a large number to be safe
    vals = np.zeros(shape=(len(etas), len(dets), 40))
    det_mats = {d: [] for d in dets}

    for i, eta in enumerate(etas):
        if i % PROGRESS_STEP == 0:
            print('Scanned {:d}/{:d} for {:.2f} <= eta <= {:.2f}'.format(i, len(etas), etas[0], etas[-1]),
                  end='\r', flush=True)
        # scan material layers
        end_x = path_r*np.cos(phi/180.*np.pi)
        end_y = path_r*np.sin(phi/180.*np.pi)
        end_z = path_r*np.sinh(eta)
        # print('({:.2f}, {:.2f}, {:.2f})'.format(end_x, end_y, end_z))
        dfr = material_scan(desc, start_point, (end_x, end_y, end_z))

        # assign material layers to detectors
        mdets = []
        for x, y, z, mat in dfr[['x', 'y', 'z', 'material']].values:
            mdet = 'Unknown'
            for det, dconf in config['detectors'].items():
                mats = dconf.get('materials', [])
                geo = dconf.get('geo', [])
                # material not match
                if not any([fnmatch.fnmatch(mat, m) for m in mats]):
                    continue
                # geo constraints
                in_geo = True
                for g in geo:
                    if not eval(g):
                        in_geo = False
                        break
                if in_geo:
                    # find the detector
                    mdet = det
                    break
            mdets.append(mdet)
        dfr.loc[:, 'detector2'] = mdets
        # print(dfr.groupby('detector')['X0'].sum())
        # aggregated values for detectors
        x0_vals = dfr.groupby('detector')['X0'].sum().to_dict()
        for j, det in enumerate(dets):
            vals[i, j, 0] = x0_vals.get(det, 0.)
            # update material dict
            dfd = dfr[dfr['detector'] == det]
            x0_mats = dfd.groupby('material')['X0'].sum()
            for mat in x0_mats.index:
                if mat not in det_mats[det]:
                    det_mats[det] = det_mats[det] + [mat]
            for k, mat in enumerate(det_mats[det]):
                vals[i, j, k + 1] = x0_mats.to_dict().get(mat, 0.)

    print('Scanned {:d}/{:d} lines for {:.2f} < eta < {:.2f}'.format(len(etas), len(etas), etas[0], etas[-1]))

    # aggregated_data
    dfa = pd.DataFrame(data=vals[:, :, 0], columns=dets, index=etas)
    # sort columns by the total thickness (over eta range)
    dets = list(dfa.sum().sort_values(ascending=False).index)
    dfa = dfa[dets]
    dfa.to_csv('material_scan_agg.csv')

    # plot
    pdf = matplotlib.backends.backend_pdf.PdfPages('material_scan_details.pdf')
    colors = ['royalblue', 'forestgreen', 'darkviolet', 'silver', 'indianred', 'goldenrod', 'darkturquoise']

    fig, ax = plt.subplots(figsize=(16, 5), dpi=160,
                           gridspec_kw={'top': 0.995, 'bottom': 0.2, 'left': 0.08, 'right': 0.98})
    bottom = np.zeros(len(dfa))
    width = np.mean(np.diff(dfa.index))
    if len(dfa.columns) > len(colors):
        print('Warning: not enough colors, ignored detector {}'.format(list(dfa.columns[len(colors):])))
    for col, c in zip(dfa.columns, colors):
        ax.fill_between(dfa.index, bottom, dfa[col].values + bottom, label=col, step='mid', color=c)
        bottom += dfa[col].values
    ax.tick_params(which='both', direction='in', labelsize=22)
    ax.set_xlabel('$\eta$', fontsize=22)
    ax.set_ylabel('X0', fontsize=22)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    ax.grid(ls=':', which='both')
    ax.set_axisbelow(False)
    ax.set_xlim(-2.2, 1.8)
    ax.set_ylim(0., ax.get_ylim()[1]*1.1)
    ax.legend(bbox_to_anchor=(0.0, 0.9, 1.0, 0.1), ncol=6, loc="upper center", fontsize=22,
          borderpad=0.2, labelspacing=0.2, columnspacing=0.6, borderaxespad=0.05, handletextpad=0.4)
    # ax.set_yscale('log')
    fig.savefig('material_scan.png')
    pdf.savefig(fig)


    for j, det in enumerate(dets):
        mats = det_mats[det]
        dfa = pd.DataFrame(data=vals[:, j, 1:len(mats)+1], columns=mats, index=etas)
        if len(dfa.columns) > len(colors):
            print('Warning: not enough colors, ignored materials {} for detector {}'.format(list(dfa.columns[len(colors):]), det))

        fig, ax = plt.subplots(figsize=(16, 6), dpi=160,
                           gridspec_kw={'top': 0.8, 'bottom': 0.2, 'left': 0.08, 'right': 0.98})
        bottom = np.zeros(len(dfa))
        width = np.mean(np.diff(dfa.index))
        for col, c in zip(dfa.columns, colors):
            ax.fill_between(dfa.index, bottom, dfa[col].values + bottom, label=col, step='mid', color=c)
            bottom += dfa[col].values
        ax.legend(bbox_to_anchor=(0.06, 1.02, 0.9, 0.1), ncol=5, loc="upper left", fontsize=18,
                  borderpad=0.3, labelspacing=0.2, columnspacing=0.8, borderaxespad=0.1, handletextpad=0.4)
        ax.tick_params(which='both', direction='in', labelsize=22)
        ax.set_xlabel('$\eta$', fontsize=22)
        ax.set_ylabel('X0', fontsize=22)
        ax.set_title(det, fontsize=22)
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.grid(ls=':', which='both')
        ax.set_axisbelow(False)
        ax.set_xlim(-2.2, 1.8)
        # ax.set_yscale('log')
        pdf.savefig(fig)
    pdf.close()
