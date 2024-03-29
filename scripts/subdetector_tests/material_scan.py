# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Chao Peng
'''
    A script to scan the materials thickness (X0 or Lambda) over eta

    It uses dd4hep::rec::MaterialManager::placementsBetween to check the material layers, and then uses the detector
    alignment to transform world coordinates to local coordiantes, and then assigns the materials to a detector based
    on TGeoVolume::Contains
    Take a grain of salt about the materials->detector assignment, especially when the detector geometry is very complex
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
OTHERS_NAME = 'Others'
# maximum number of detectors/materials to plot
COLORS = ['royalblue', 'indianred', 'forestgreen', 'gold', 'darkviolet', 'orange', 'darkturquoise']
# a specified color for Others, this should not be included in COLORS
OTHERS_COLOR = 'silver'


# FIXME: this is a work-around to an issue of dd4hep::rec::MaterialManager::placementsBetween
# As of 04/21/2023, at negative eta (eta < 0.), the scan will miss the first material layer (vacuum, 2.8 cm)
# To reporduce this issue, check the difference between:
#   materialScan epic_craterlake.xml 0 0 0 100  40  -0.01 | grep Vacuum
#   materialScan epic_craterlake.xml 0 0 0 100  40  0.01 | grep Vacuum
# It is GEOMETRY DEPENDENT (the thickness of the first vacuum material layer)
# This helper class can be removed once the issue is fixed
# Plan to file a bug report for this (as of 04/21/2023)
class ThicknessCorrector:
    def __init__(self, desc=None):
        self.missing = 0.
        self.scanned = False
        if desc is not None:
            self.scan_missing_thickness(desc)

    def scan_missing_thickness(self, desc, pi=(0.,0.,0.), pf=(100.,40.,0.), dz=0.01, epsilon=1e-3):
        # assume negative eta is missed, maybe can do a more detailed check?
        pf1 = np.array(pf) + np.array([0., 0., dz])
        dft1 = material_scan(desc, pi, pf1, epsilon)
        th1 = dft1['path_length'].iloc[0]
        pf2 = np.array(pf) + np.array([0., 0., -dz])
        dft2 = material_scan(desc, pi, pf2, epsilon)
        th2 = dft2['path_length'].iloc[0]
        self.scanned = True
        # print(dft1.head(3))
        # print(dft2.head(3))
        if th1 != th2:
            self.missing = th1

    def correct_path_length(self, direction, pl=0.):
        if direction[2] < 0.:
            pl2 = pl + np.sqrt(1./(direction[0]**2 + direction[1]**2))*self.missing
            # print('correcting path length {:.2f} -> {:.2f}'.format(pl, pl2))
            return pl2
        return pl



'''
    Re-implementation of MaterialScan::Print in DD4Hep::DDRec
    MaterialScan does not have a python interface and that function is relatively simple, so here it is
    desc: DD4hep::Detector
    start: 3D vector for start point
    end: 3D vector for end point
    epsilon: step size
'''
def material_scan(desc, start, end, epsilon, int_dets=None, thickness_corrector=None):
    mat_mng = DDRec.MaterialManager(desc.worldVolume())
    # only use the top-level detectors
    if int_dets is None:
        dets_list = [d for n, d in desc.world().children()]
    else:
        dets_list = [d for n, d in desc.world().children() if n in int_dets]
    # rvec = ROOT.Math.XYZVector()
    # id_conv = DDRec.CellIDPositionConverter(desc)
    # det_dict = {d.id(): n for n, d in desc.world().children()}

    p0 = np.array(start)
    p1 = np.array(end)
    direction = (p1 - p0)/np.linalg.norm(p1 - p0)
    # print(p0, p1, direction)
    placements = mat_mng.placementsBetween(tuple(p0), tuple(p1), epsilon);

    # calculate material layer by layer
    int_x0 = 0
    int_lambda = 0
    path_length = 0.
    if thickness_corrector is not None:
        path_length = thickness_corrector.correct_path_length(direction, path_length)

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
        d0 = OTHERS_NAME
        for d in dets_list:
            local = d.nominal().worldToLocal(pcurr)
            local = np.array([local.X(), local.Y(), local.Z()])
            if d.volume().Contains(local):
                d0 = d.GetName()
        res.append([
            d0,
            # det_dict.get(det_id, 'Unknown'),
            mat.GetName(), mat.GetZ(), mat.GetA(), mat.GetDensity(),
            radl, intl, l, path_length,
            x0, lmd,
            pcurr[0], pcurr[1], pcurr[2],
            # local[0], local[1], local[2],
            ])
    cols = [
        'detector',
        'material', 'Z', 'A', 'density',
        'radl', 'intl', 'thickness', 'path_length',
        'X0', 'lambda',
        'x', 'y', 'z',
        # 'local_x', 'local_y', 'local_z'
        ]
    dft = pd.DataFrame(data=res, columns=cols)
    # print(dft.groupby('detector')['X0'].sum())
    return dft


# the allowed value types {name: [col_name, label_name, major_step, minor_step]}
ALLOWED_VALUE_TYPES = {
    'X0': ['X0', '$X_0$', 10, 5],
    'lambda': ['lambda', '$\lambda$', 5, 1],
}
# execute the script
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            dest='compact',
            help='Top-level xml file of the detector description.'
            )
    parser.add_argument(
            '--path-r', type=float, default=400., # 120.,
            help='R_xy (cm) where the scan stops.'
            )
    parser.add_argument(
            '--start-point', default='0,0,0',
            help='Start point of the scan, use the format \"x,y,z\", unit is cm.'
            )
    parser.add_argument(
            '--eta-min', type=float, default=-2.0,
            help='Minimum eta for the scan.'
            )
    parser.add_argument(
            '--eta-max', type=float, default=2.0,
            help='Minimum eta for the scan.'
            )
    parser.add_argument(
            '--eta-nbins', type=int, default=401,
            help='Number of bins for the eta scan.'
            )
    parser.add_argument(
            '--phi', type=float, default=20.,
            help='Phi angle of the scan, unit is degree.'
            )
    parser.add_argument(
            '--epsilon', type=float, default=1e-3,
            help='Step size for each material scan.'
            )
    parser.add_argument(
            '--value-type', default='X0',
            help='Choose one in {}.'.format(list(ALLOWED_VALUE_TYPES.keys()))
            )
    parser.add_argument(
            '--detectors',
            # default='all',
            # ordering the detectors
            default='EcalBarrelScFi,EcalEndcapN,EcalEndcapP,SolenoidBarrel,HcalBarrel,HcalEndcapN,HcalEndcapP',
            help='Names of the interested detectors, separated by \",\". Use \"all\" for all detectors.'
            )
    args = parser.parse_args()

    if not os.path.exists(args.compact):
        print('Cannot find {}'.format(args.compact))
        exit(-1)

    if args.value_type not in ALLOWED_VALUE_TYPES.keys():
        print('Cannot find value {}, choose one in {}'.format(args.value_type, list(ALLOWED_VALUE_TYPES.keys())))
        exit(-1)

    # scan parameters
    path_r = args.path_r
    phi = args.phi
    etas = np.linspace(args.eta_min, args.eta_max, args.eta_nbins)
    start_point = np.array([float(v.strip()) for v in args.start_point.split(',')])
    if len(start_point) != 3:
        print('Error: expecting three values (x,y,z) from --start-point, getting {:d} instead.'.format(len(start_point)))
        exit(-1)

    # geometry initialization
    desc = dd4hep.Detector.getInstance()
    desc.fromXML(args.compact)

    # check detector list (case sensitive)
    all_dets = [n for n, _ in desc.world().children()]
    dets = []
    missing_dets = []
    sort_dets = False
    for det in args.detectors.split(','):
        det = det.strip()
        if det.lower() == 'all':
            dets = all_dets
            # sort detectors by material thickness
            sort_dets = True
            break
        if det in all_dets:
            dets.append(det)
        else:
            missing_dets.append(det)
    if len(missing_dets) > 0:
        print('Warning: detectors {} were missing from the description.'.format(missing_dets))
        print('A full list of detectors are:')
        for d in all_dets:
            print('   --- {}'.format(d))


    # FIXME: work-around for dd4hep material scan issue, check ThicknessCorrector
    th_corr = ThicknessCorrector()
    th_corr.scan_missing_thickness(desc, start_point, np.array(start_point) + np.array([100., 50., 0.]),
                                   dz=0.0001, epsilon=args.epsilon)

    # array for detailed data
    # number of materials cannot be pre-determined, so just assign a large number to be safe
    vals = np.zeros(shape=(len(etas), len(dets) + 1, 50))
    dets2 = dets + [OTHERS_NAME]
    det_mats = {d: [] for d in dets2}
    vt = ALLOWED_VALUE_TYPES[args.value_type]

    for i, eta in enumerate(etas):
        if i % PROGRESS_STEP == 0:
            print('Scanned {:d}/{:d} for {:.2f} <= eta <= {:.2f}'.format(i, len(etas), etas[0], etas[-1]),
                  end='\r', flush=True)
        # scan material layers
        end_x = path_r*np.cos(phi/180.*np.pi)
        end_y = path_r*np.sin(phi/180.*np.pi)
        end_z = path_r*np.sinh(eta)
        # print('({:.2f}, {:.2f}, {:.2f})'.format(end_x, end_y, end_z))
        dfr = material_scan(desc, start_point, (end_x, end_y, end_z), epsilon=args.epsilon, int_dets=dets, thickness_corrector=th_corr)

        # aggregated values for detectors
        x0_vals = dfr.groupby('detector')[vt[0]].sum().to_dict()
        for j, det in enumerate(dets2):
            vals[i, j, 0] = x0_vals.get(det, 0.)
            # update material dict
            dfd = dfr[dfr['detector'] == det]
            x0_mats = dfd.groupby('material')[vt[0]].sum()
            for mat in x0_mats.index:
                if mat not in det_mats[det]:
                    det_mats[det] = det_mats[det] + [mat]
            for k, mat in enumerate(det_mats[det]):
                vals[i, j, k + 1] = x0_mats.to_dict().get(mat, 0.)

    print('Scanned {:d}/{:d} lines for {:.2f} < eta < {:.2f}'.format(len(etas), len(etas), etas[0], etas[-1]))

    # aggregated data for plots
    dfa = pd.DataFrame(data=vals[:, :, 0], columns=dets2, index=etas)
    # sort columns by the total thickness (over eta range)
    sdets = dets
    if sort_dets:
        sdets = [d for d in dfa.sum().sort_values(ascending=False).index if d != OTHERS_NAME]
        dfa = dfa[sdets + [OTHERS_NAME]]
    dfa.to_csv('material_scan_agg.csv')

    # plots
    pdf = matplotlib.backends.backend_pdf.PdfPages('material_scan_details.pdf')
    fig, ax = plt.subplots(figsize=(16, 8), dpi=160,
                           gridspec_kw={'top': 0.9, 'bottom': 0.2, 'left': 0.08, 'right': 0.98})

    # group some detectors into Others because not enough colors for them
    plot_dets = sdets
    if len(sdets) > len(COLORS):
        group_dets = sdets[len(COLORS):]
        dfa.loc[:, OTHERS_NAME] += dfa.loc[:, group_dets].sum(axis=1)
        print('Detectors {} are grouped into {} due to limited number of colors.'.format(group_dets, OTHERS_NAME))
        plot_dets = sdets[:len(COLORS)]

    # plot
    bottom = np.zeros(len(dfa))
    width = np.mean(np.diff(dfa.index))
    for col, c in zip(plot_dets, COLORS):
        ax.fill_between(dfa.index, bottom, dfa[col].values + bottom, label=col, step='mid', color=c)
        bottom += dfa[col].values
    # only plot Others if there are any
    if dfa[OTHERS_NAME].sum() > 0.:
        ax.fill_between(dfa.index, bottom, dfa[OTHERS_NAME].values + bottom, label=OTHERS_NAME, step='mid', color=OTHERS_COLOR)

    # formatting
    ax.tick_params(which='both', direction='in', labelsize=22)
    ax.set_xlabel('$\eta$', fontsize=22)
    ax.set_ylabel('{} ($R_{{xy}} \leq {:.1f}$ cm)'.format(vt[1], args.path_r), fontsize=22)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(vt[2]))
    ax.yaxis.set_minor_locator(MultipleLocator(vt[3]))
    ax.grid(ls=':', which='both')
    ax.set_axisbelow(False)
    ax.set_xlim(args.eta_min, args.eta_max)
    ax.set_ylim(0., ax.get_ylim()[1]*1.1)
    ax.legend(bbox_to_anchor=(0.0, 0.9, 1.0, 0.1), ncol=4, loc="upper center", fontsize=22,
          borderpad=0.2, labelspacing=0.2, columnspacing=0.6, borderaxespad=0.05, handletextpad=0.4)
    # ax.set_yscale('log')
    fig.savefig('material_scan.png')
    pdf.savefig(fig)
    plt.close(fig)

    # detailed plot for each detector
    for j, det in enumerate(dets):
        dmats = det_mats[det]
        dfa = pd.DataFrame(data=vals[:, j, 1:len(dmats)+1], columns=dmats, index=etas)
        if not len(dmats):
            print('No material found for detector {} in this scan, skipped it.'.format(det))
            continue
        # sort columns by the total thickness (over eta range)
        mats = [d for d in dfa.sum().sort_values(ascending=False).index]
        dfa = dfa[mats]
        plot_mats = mats
        if len(mats) > len(COLORS):
            group_mats = mats[len(COLORS):]
            dfa.loc[:, OTHERS_NAME] = dfa.loc[:, group_mats].sum(axis=1)
            print('For detector {}, materials {} are grouped into {} due to limited number of colors.'.format(det, group_mats, OTHERS_NAME))
            plot_mats = mats[:len(COLORS)]

        fig, ax = plt.subplots(figsize=(16, 8), dpi=160,
                           gridspec_kw={'top': 0.9, 'bottom': 0.2, 'left': 0.08, 'right': 0.98})
        bottom = np.zeros(len(dfa))
        width = np.mean(np.diff(dfa.index))
        for col, c in zip(plot_mats, COLORS):
            ax.fill_between(dfa.index, bottom, dfa[col].values + bottom, label=col, step='mid', color=c)
            bottom += dfa[col].values
        if OTHERS_NAME in dfa.columns and dfa[OTHERS_NAME].sum() > 0.:
            ax.fill_between(dfa.index, bottom, dfa[OTHERS_NAME].values + bottom, label=OTHERS_NAME, step='mid', color=OTHERS_COLOR)

        ax.legend(bbox_to_anchor=(0.0, 0.9, 1.0, 0.1), ncol=6, loc="upper left", fontsize=18,
                  borderpad=0.3, labelspacing=0.2, columnspacing=0.8, borderaxespad=0.1, handletextpad=0.4)
        ax.tick_params(which='both', direction='in', labelsize=22)
        ax.set_xlabel('$\eta$', fontsize=22)
        ax.set_ylabel('{} ($R_{{xy}} \leq {:.1f}$ cm)'.format(vt[1], args.path_r), fontsize=22)
        ax.set_title(det, fontsize=22)
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(vt[2]))
        ax.yaxis.set_minor_locator(MultipleLocator(vt[3]))
        ax.grid(ls=':', which='both')
        ax.set_axisbelow(False)
        ax.set_xlim(args.eta_min, args.eta_max)
        ax.set_ylim(0., ax.get_ylim()[1]*1.1)
        pdf.savefig(fig)
        plt.close(fig)
    # close PDF
    pdf.close()
