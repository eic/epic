#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 Sylvester Joosten
'''
    Draw x-radial cross-sections of the BIC imaging layers (staves, carbon plates and side walls)
    from the loaded DD4hep/TGeo geometry. x runs across the sector, z radially outward from the
    inner face of the layer; layer 1 is the bottom panel of the combined figure.
    Staves are blue (even index, low rail) or orange (odd, high rail), the sensitive Si slice red,
    carbon plates and walls grey, the layer trapezoid black.

    use case (writes <prefix>_layer<N>.png, <prefix>_all.png and <prefix>_zoom.png):
    python scripts/subdetector_tests/bic/draw_bic_imaging_sector.py -c epic_craterlake_bic_6layers.xml \
        -o figures/after --label After --layers 1-6
'''
import argparse
import array
import math
import os

import ROOT
import dd4hep
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

ROOT.gROOT.SetBatch(True)

# volume-name prefixes defined by src/BarrelCalorimeterImaging_geo.cpp
SECTOR = 'sector'
LAYER = 'layer'
STAVE = 'stave'
PLATE = 'frame'            # inner/outer carbon plates (frame_inner<N>, frame_outer<N>)
WALL = 'frame_side'        # carbon side walls, checked before PLATE

LOW, HIGH, SI, GREY = '#2166ac', '#e08214', '#b2182b', '#777777'


def box_pts(sh):
    return [(-sh.GetDX(), -sh.GetDZ()), (sh.GetDX(), -sh.GetDZ()),
            (sh.GetDX(), sh.GetDZ()), (-sh.GetDX(), sh.GetDZ())]


def trd_pts(sh):
    return [(-sh.GetDx1(), -sh.GetDz()), (sh.GetDx1(), -sh.GetDz()),
            (sh.GetDx2(), sh.GetDz()), (-sh.GetDx2(), sh.GetDz())]


def trap_pts(sh):
    # TGeoTrap is a TGeoArb8: 4 (x, y) vertices at -dz, then 4 at +dz; project on x-z
    v = sh.GetVertices()
    dz = sh.GetDz()
    lo = [v[2 * k] for k in range(4)]
    hi = [v[2 * k] for k in range(4, 8)]
    return [(min(lo), -dz), (max(lo), -dz), (max(hi), dz), (min(hi), dz)]


def shape_pts(sh):
    return trap_pts(sh) if sh.IsA().GetName() == 'TGeoTrap' else trd_pts(sh)


def xf(mats, pts, zoff):
    '''transform shape points through the node matrices (innermost first), cm -> mm, z from the layer bottom'''
    out = []
    for x, z in pts:
        p = array.array('d', [x, 0., z])
        for m in mats:
            q = array.array('d', [0., 0., 0.])
            m.LocalToMaster(p, q)
            p = q
        out.append((10 * p[0], 10 * (p[2] + zoff)))
    return out


def find_sectors(vol, names, depth=0):
    '''sector volumes below the detector volume, by name'''
    found = {}
    for i in range(vol.GetNdaughters()):
        dv = vol.GetNode(i).GetVolume()
        if dv.GetName() in names:
            found[dv.GetName()] = dv
        elif depth < 4:
            found.update(find_sectors(dv, names, depth + 1))
    return found


def extract(args):
    det = dd4hep.Detector.getInstance()
    det.fromXML(args.compact)
    sdet = det.detector(args.detector)
    if not sdet.isValid():
        raise SystemExit('detector %s not found in %s' % (args.detector, args.compact))
    secs = find_sectors(sdet.placement().volume().ptr(), set(args.sector))
    if not secs:
        raise SystemExit('none of the volumes %s found in %s' % (args.sector, args.detector))
    layers = {}
    for sv in secs.values():
        for i in range(sv.GetNdaughters()):
            ln = sv.GetNode(i)
            if not ln.GetName().startswith(LAYER):
                continue
            n = int(ln.GetName()[len(LAYER):].split('_')[0])
            if n not in args.layers:
                continue
            lv = ln.GetVolume()
            trd = lv.GetShape()
            zoff = trd.GetDz()  # layer-local z measured from the inner face
            d = dict(outline=xf([], trd_pts(trd), zoff), plates=[], walls=[], staves=[])
            for j in range(lv.GetNdaughters()):
                dn = lv.GetNode(j)
                name = dn.GetName()
                sh = dn.GetVolume().GetShape()
                if name.startswith(WALL):
                    d['walls'].append(xf([dn.GetMatrix()], shape_pts(sh), zoff))
                elif name.startswith(PLATE):
                    d['plates'].append(xf([dn.GetMatrix()], shape_pts(sh), zoff))
                elif name.startswith(STAVE):
                    st = dn.GetVolume()
                    mn = st.GetNode(0)
                    cn = mn.GetVolume().GetNode(0)
                    sn = cn.GetVolume().GetNode(0)  # sensitive slice
                    d['staves'].append(dict(
                        idx=len(d['staves']),  # placement order = stave index within the layer
                        box=xf([dn.GetMatrix()], box_pts(sh), zoff),
                        si=xf([sn.GetMatrix(), cn.GetMatrix(), mn.GetMatrix(), dn.GetMatrix()],
                              box_pts(sn.GetVolume().GetShape()), zoff)))
            if d['staves']:
                layers[n] = d
    return layers


def wall_slope(d):
    '''half-angle of the slanted layer side wall: atan(dx/dz) of the trapezoid outline'''
    (xa, za), (xb, zb) = d['outline'][0], d['outline'][3]
    return math.atan2(abs(xb - xa), zb - za)


def draw_layer(ax, d, n, label, stretch=True):
    ax.add_patch(Polygon(d['outline'], closed=True, fill=False, ec='k', lw=1.2))
    for pl in d['plates']:
        ax.add_patch(Polygon(pl, closed=True, fc=GREY, ec=GREY, alpha=0.6))
    for pl in d['walls']:
        ax.add_patch(Polygon(pl, closed=True, fc='#4d4d4d', ec='#4d4d4d', alpha=0.8))
    for s in d['staves']:
        col = LOW if s['idx'] % 2 == 0 else HIGH
        ax.add_patch(Polygon(s['box'], closed=True, fc=col, ec=col, alpha=0.35, lw=0.8))
        ax.add_patch(Polygon(s['si'], closed=True, fc=SI, ec=SI, alpha=0.9))
    xs = [x for x, _ in d['outline']]
    zs = [z for _, z in d['outline']]
    ax.set_xlim(min(xs) - 3, max(xs) + 3)
    ax.set_ylim(min(zs) - 1, max(zs) + 1)
    st = d['staves']
    pitch = st[1]['box'][0][0] - st[0]['box'][0][0] if len(st) > 1 else 0
    w = st[0]['box'][1][0] - st[0]['box'][0][0]
    ax.set_title('%slayer %d: %d staves, pitch %.2f mm, stave overlap %.2f mm'
                 % (label + ', ' if label else '', n, len(st), pitch, w - pitch), fontsize=10)
    ax.set_ylabel('radial z [mm] (up = outward)')
    ax.set_xlabel('x across the sector [mm]' + ('  (vertical axis stretched)' if stretch else ''))


def zoom(d, n, label, out):
    '''three adjacent staves at the left wall, true scale, annotated'''
    st = d['staves']
    fig, ax = plt.subplots(figsize=(11, 6.5))
    draw_layer(ax, d, n, label, stretch=False)
    x0, zb = d['outline'][0]
    s0, s1, s2 = (s['box'] for s in st[:3])
    xc0 = 0.5 * (s0[0][0] + s0[1][0])
    xc1 = 0.5 * (s1[0][0] + s1[1][0])
    top = max(p[1] for b in (s0, s1, s2) for p in b) - zb  # height of the stave stack above the layer bottom
    ax.set_xlim(x0 - 2, s1[1][0] + 0.6 * (xc1 - xc0))  # up to just past the second stave
    ax.set_ylim(zb - 1, zb + top + 5)
    ax.set_aspect('equal')
    ax.set_title('%slayer %d, three adjacent staves at the sector wall (true scale)'
                 % (label + ': ' if label else '', n), fontsize=11)
    tilted = abs(s0[0][1] - s0[1][1]) > 1e-6
    ax.annotate('', xy=(xc1, zb + top + 2), xytext=(xc0, zb + top + 2), arrowprops=dict(arrowstyle='<->'))
    ax.text(0.5 * (xc0 + xc1), zb + top + 2.4, 'pitch %.2f mm' % (xc1 - xc0), ha='center', fontsize=9)
    if not tilted:
        ov0, ov1 = s1[0][0], s0[1][0]
        ax.annotate('', xy=(ov1, zb + top + 0.4), xytext=(ov0, zb + top + 0.4),
                    arrowprops=dict(arrowstyle='<->', color='r'))
        ax.text(0.5 * (ov0 + ov1), zb + top + 0.8, 'overlap %.2f mm' % (ov1 - ov0), ha='center', fontsize=9, color='r')
        # rails: heights of stave bottoms above the plate outer face
        plate_top = max(p[1] for p in d['plates'][0])
        for s in st[:2]:
            b = s['box']
            xr = b[0][0] + 1.5
            ax.annotate('', xy=(xr, b[0][1]), xytext=(xr, plate_top), arrowprops=dict(arrowstyle='<->'))
            ax.text(xr + 0.4, 0.5 * (b[0][1] + plate_top), 'rail %.2f mm' % (b[0][1] - plate_top),
                    fontsize=8, va='center')
        if d['walls']:  # inner face of the left wall: the max-x edge of the wall polygon
            w = min(d['walls'], key=lambda q: q[0][0])
            (xa, za), (xb, zc) = w[1], w[2]

            def wall_x(z):
                return xa + (z - za) * (xb - xa) / (zc - za)
            what = 'wall inner face'
        else:
            (xa, za), (xb, zc) = d['outline'][0], d['outline'][3]

            def wall_x(z):
                return xa + (z - za) * (xb - xa) / (zc - za)
            what = 'wall'
        c = (s0[0][0] - wall_x(s0[0][1])) * math.cos(wall_slope(d))  # perpendicular to the slanted face
        ax.annotate('min. distance to %s\n(perpendicular) %.2f mm' % (what, abs(c)), xy=s0[0],
                    xytext=(x0 + 1, zb + top - 0.5), fontsize=8, arrowprops=dict(arrowstyle='->'))
    ax.text(x0 + 1, zb + top + 3.8,
            'blue: even stave index%s; orange: odd stave index%s; red: sensitive Si; grey: carbon plates and side walls'
            % ('' if tilted else ' (low rail)', '' if tilted else ' (high rail)'), fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def parse_layers(text):
    '''"1,3-5" -> [1, 3, 4, 5]'''
    out = []
    for part in text.split(','):
        lo, _, hi = part.partition('-')
        out += range(int(lo), int(hi or lo) + 1)
    return sorted(set(out))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='draw_bic_imaging_sector.py',
        description='Draw x-radial cross-sections of the BIC imaging layers from the loaded geometry')
    parser.add_argument('-c', '--compact', help='compact detector file',
                        default=f'{os.getenv("DETECTOR_PATH")}/{os.getenv("DETECTOR_CONFIG")}.xml')
    parser.add_argument('-o', '--output', default='bic_imaging', help='output file prefix')
    parser.add_argument('--detector', default='EcalBarrelImaging', help='name of the imaging calorimeter')
    parser.add_argument('--sector', nargs='+', default=['sector1', 'sector26'],
                        help='sector volumes to take layers from (layer 1 lives in sector1, layers 2-6 in sector26)')
    parser.add_argument('--layers', default='1', help='layers to draw, e.g. 1 or 1-6 or 1,3')
    parser.add_argument('--zoom-layer', type=int, help='layer for the zoom figure (default: first drawn layer)')
    parser.add_argument('--label', default='', help='label for the figure titles')
    args = parser.parse_args()
    args.layers = set(parse_layers(args.layers))
    if args.zoom_layer:
        args.layers.add(args.zoom_layer)

    layers = extract(args)
    if not layers:
        raise SystemExit('no layers %s with staves found' % sorted(args.layers))
    for n in sorted(layers):
        fig, ax = plt.subplots(figsize=(14, 3.4))
        draw_layer(ax, layers[n], n, args.label)
        fig.tight_layout()
        fig.savefig('%s_layer%d.png' % (args.output, n), dpi=130)
        plt.close(fig)
        print('wrote %s_layer%d.png' % (args.output, n))
    if len(layers) > 1:
        fig, axs = plt.subplots(len(layers), 1, figsize=(14, 2.9 * len(layers)))
        for ax, n in zip(axs, sorted(layers, reverse=True)):  # layer 1 at the bottom
            draw_layer(ax, layers[n], n, args.label)
        fig.tight_layout()
        fig.savefig(args.output + '_all.png', dpi=110)
        plt.close(fig)
        print('wrote %s_all.png' % args.output)
    zl = args.zoom_layer or min(layers)
    if zl in layers and len(layers[zl]['staves']) >= 3:
        zoom(layers[zl], zl, args.label, args.output + '_zoom.png')
        print('wrote %s_zoom.png' % args.output)
    else:
        print('no zoom: layer %d not found or has fewer than 3 staves' % zl)
