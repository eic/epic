#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 Sylvester Joosten
'''
    3D render (orthographic z-buffer) of a few BIC imaging layers from the loaded DD4hep/TGeo geometry.
    Walks the TGeo nodes of a sector volume, clips them to a window along the beam (y) and radially (z,
    sector frame), skips the outer carbon plates and the ScFi/Pb volumes (--skip, a safety net when a sector holding them is selected), and draws every leaf box/trd/trap
    as shaded faces: stave tops in the stave tier colour (blue even index, orange odd), red Si, yellow
    Kapton, orange Cu underneath, dark grey carbon plates and side walls.

    use case (layers 2 and 3 of sector26, +-3 cm along the beam, radial direction up):
    python scripts/subdetector_tests/bic/render_bic_imaging_3d.py -c epic_craterlake_bic_6layers.xml \
        -o render3d_layers2-3.png --sector sector26 --y-range -3 3 --view 24 62 --z-range 86.7 92.65 \
        --title 'BIC imaging layers 2 (bottom) and 3 (top), radial direction up; outer carbon plates and ScFi/Pb hidden'
'''
import argparse
import array
import os
import re

import ROOT
import dd4hep
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

ROOT.gROOT.SetBatch(True)

# volume-name prefixes defined by src/BarrelCalorimeterImaging_geo.cpp
LAYER = 'layer'
STAVE = 'stave'
OUTER_PLATE = 'frame_outer'

SS = 3  # supersampling
COLS = {'Silicon': '#c0392b', 'CarbonFiber': '#555555', 'Copper': '#e67e22', 'Kapton': '#f1c40f', 'Epoxy': '#27ae60'}
TIER = ('#2166ac', '#e08214')  # even stave index: low rail, odd: high rail (as in draw_bic_imaging_sector.py)
FACES = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]


def shape_pts(shape):
    a = array.array('d', [0.] * 24)
    shape.SetPoints(a)
    return np.array(a).reshape(8, 3)


def to_master(mats, p):
    for mm in reversed(mats):
        b = array.array('d', [0.] * 3)
        mm.LocalToMaster(array.array('d', p), b)
        p = b
    return np.array(p)


def collect(vol, args, skip):
    '''leaf box/trd/trap volumes in the window, as (8 corner points, colour) pairs'''
    ylo, yhi = args.y_range
    zlo, zhi = args.z_range
    polys, fc, count = [], [], {}
    tier = [0]
    stave_idx = [-1]

    def walk(v, mats, depth):
        for i in range(v.GetNdaughters()):
            n = v.GetNode(i)
            name = n.GetName()
            dv = n.GetVolume()
            mat = dv.GetMaterial().GetName()
            if skip.search(name) or skip.search(mat) or name.startswith(OUTER_PLATE):
                continue
            m = mats + [n.GetMatrix()]
            if depth == 0 and name.startswith(LAYER):
                stave_idx[0] = -1
            if depth == 1 and name.startswith(STAVE):
                stave_idx[0] += 1
                tier[0] = stave_idx[0] % 2  # even: low rail
            o = to_master(m, [0., 0., 0.])
            if not (zlo - 3 < o[2] < zhi + 3 or depth < 1):
                continue
            if depth >= 2 and (o[1] > yhi + 12 or o[1] < ylo - 12):
                continue
            sh = dv.GetShape()
            if mat in COLS and dv.GetNdaughters() == 0 and (
                    sh.IsA().InheritsFrom('TGeoBBox') or sh.IsA().GetName() in ('TGeoTrd2', 'TGeoTrap')):
                # leaf volumes only: the chip container is Si too
                P = np.array([to_master(m, list(q)) for q in shape_pts(sh)])
                if zlo < P[:, 2].mean() < zhi:
                    P[:, 1] = np.clip(P[:, 1], ylo, yhi)
                    if P[:, 1].max() - P[:, 1].min() > 1e-6:
                        c = TIER[tier[0]] if mat == 'Epoxy' else COLS[mat]  # stave top: tier colour
                        count[mat] = count.get(mat, 0) + 1
                        for f in FACES:
                            polys.append(P[f])
                            fc.append(matplotlib.colors.to_rgb(c))
            walk(dv, m, depth + 1)
    walk(vol, [], 0)
    print('boxes:', count, 'faces', len(polys))
    return polys, fc


def render(polys, fc, elev, azim):
    # orthographic camera, radial z up on screen
    el, az = np.radians(elev), np.radians(azim)
    D = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])  # towards the viewer
    Rt = np.array([-np.sin(az), np.cos(az), 0.])
    Up = np.cross(D, Rt)
    allp = np.vstack(polys)
    sx, sy = allp @ Rt, allp @ Up
    W, H = 1800 * SS, 1100 * SS
    sc = min(0.94 * W / np.ptp(sx), 0.94 * H / np.ptp(sy))
    ox, oy = W / 2 - sc * (sx.max() + sx.min()) / 2, H / 2 + sc * (sy.max() + sy.min()) / 2
    zbuf = np.full((H, W), -np.inf)
    img = np.ones((H, W, 3))
    oid = np.full((H, W), -1, dtype=np.int32)
    light = np.array([0.35, -0.45, 0.82])
    light /= np.linalg.norm(light)

    def raster(tri, depth, colour, ident):
        xs, ys = tri[:, 0], tri[:, 1]
        x0, x1 = int(max(xs.min(), 0)), int(min(xs.max() + 1, W - 1))
        y0, y1 = int(max(ys.min(), 0)), int(min(ys.max() + 1, H - 1))
        if x1 <= x0 or y1 <= y0:
            return
        X, Y = np.meshgrid(np.arange(x0, x1) + 0.5, np.arange(y0, y1) + 0.5)
        d = (ys[1] - ys[2]) * (xs[0] - xs[2]) + (xs[2] - xs[1]) * (ys[0] - ys[2])
        if abs(d) < 1e-12:
            return
        a = ((ys[1] - ys[2]) * (X - xs[2]) + (xs[2] - xs[1]) * (Y - ys[2])) / d
        b = ((ys[2] - ys[0]) * (X - xs[2]) + (xs[0] - xs[2]) * (Y - ys[2])) / d
        c = 1 - a - b
        m = (a >= -1e-9) & (b >= -1e-9) & (c >= -1e-9)
        z = a * depth[0] + b * depth[1] + c * depth[2]
        sub = zbuf[y0:y1, x0:x1]
        upd = m & (z > sub)
        sub[upd] = z[upd]
        img[y0:y1, x0:x1][upd] = colour
        oid[y0:y1, x0:x1][upd] = ident

    for k, (P, col) in enumerate(zip(polys, fc)):
        n = np.cross(P[1] - P[0], P[2] - P[1])
        n /= (np.linalg.norm(n) + 1e-30)
        if np.dot(n, D) < 0:
            n = -n  # all solids: shade the visible side
        shade = 0.45 + 0.55 * max(np.dot(n, light), 0.) * 0.9 + 0.08 * max(np.dot(n, D), 0.)
        colour = np.clip(np.array(col) * shade, 0, 1)
        scr = np.stack([ox + sc * (P @ Rt), oy - sc * (P @ Up)], axis=1)
        dep = P @ D
        for tri in ((0, 1, 2), (0, 2, 3)):
            raster(scr[list(tri)], dep[list(tri)], colour, k)
    # outline where the face id changes
    edge = np.zeros((H, W), bool)
    edge[:, 1:] |= oid[:, 1:] != oid[:, :-1]
    edge[1:, :] |= oid[1:, :] != oid[:-1, :]
    img[edge] = img[edge] * 0.55
    return img.reshape(H // SS, SS, W // SS, SS, 3).mean(axis=(1, 3))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='render_bic_imaging_3d.py',
        description='3D z-buffer render of BIC imaging layers from the loaded geometry')
    parser.add_argument('-c', '--compact', help='compact detector file',
                        default=f'{os.getenv("DETECTOR_PATH")}/{os.getenv("DETECTOR_CONFIG")}.xml')
    parser.add_argument('-o', '--output', default='bic_imaging_3d.png', help='output figure')
    parser.add_argument('--sector', default='sector26', help='sector volume to draw (layers 2-6: sector26)')
    parser.add_argument('--y-range', nargs=2, type=float, default=[-3., 3.], metavar=('LO', 'HI'),
                        help='clip window along the beam [cm]')
    parser.add_argument('--z-range', nargs=2, type=float, required=True, metavar=('LO', 'HI'),
                        help='radial window in the sector frame [cm], selects the layers')
    parser.add_argument('--view', nargs=2, type=float, default=[24., 62.], metavar=('ELEV', 'AZIM'),
                        help='camera elevation and azimuth [deg]')
    parser.add_argument('--skip', default='scifi|^fiber',
                        help='regex (case-insensitive) on volume and material names to skip with their daughters')
    parser.add_argument('--title', default='BIC imaging layers, radial direction up; outer carbon plates and ScFi/Pb hidden')
    args = parser.parse_args()

    det = dd4hep.Detector.getInstance()
    det.fromXML(args.compact)
    vols = [v for v in ROOT.gGeoManager.GetListOfVolumes() if v.GetName() == args.sector]
    if not vols:
        raise SystemExit('volume %s not found in %s' % (args.sector, args.compact))
    polys, fc = collect(vols[0], args, re.compile(args.skip, re.I))
    if not polys:
        raise SystemExit('nothing to draw in the given windows')
    im = render(polys, fc, *args.view)

    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_axes([0, 0.05, 1, 0.9])
    ax.imshow(im)
    ax.axis('off')
    leg = [Patch(fc=c, label=l) for l, c in (
        ('stave top, low rail (even index)', TIER[0]), ('stave top, high rail (odd index)', TIER[1]),
        ('Si sensor, Cu, Kapton underneath (facing the inner plate)', COLS['Silicon']),
        ('carbon: inner plate and side walls', COLS['CarbonFiber']))]
    ax.legend(handles=leg, loc='upper left', fontsize=11, frameon=False)
    fig.suptitle(args.title, fontsize=12)
    fig.savefig(args.output, dpi=110)
    print('wrote', args.output)
