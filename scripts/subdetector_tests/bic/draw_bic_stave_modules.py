#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 Sylvester Joosten
'''
    Draw the silicon layout of one BIC imaging stave position (staveN and staveP, 12 AstroPix modules
    each, 9 chips per module) from the loaded DD4hep/TGeo geometry, seen from above: y runs along the
    beam (stave length, true scale), x across the stave (not to scale in the full-stave panel).
    All dimensions are measured from the node shapes and matrices; the gaps are printed to stdout.

    Figure 1 (<prefix>.png): the full stave position plus four equal-scale zooms (outer stave end,
    module-to-module gap, chip-to-chip gap, staveN-staveP junction at the calorimeter centre).
    Figure 2 (<prefix>_module.png): module 0 of staveN and of staveP with its 9 chips at true scale.
    Module and chip labels are the readout ids (physVolIDs) of the placements; positions go through the
    full node matrices (staveP is rotated by 180 degrees), ids count from the outer stave end.
    Colour = module id parity (blue even, orange odd); the same colours are used in all panels.

    use case (writes <prefix>.png and <prefix>_module.png):
    python scripts/subdetector_tests/bic/draw_bic_stave_modules.py -c epic_craterlake_bic_6layers.xml \
        -o pr-material/figures/bic_stave_modules --sector sector26 --layer 2 --position 0
'''
import argparse
import array
import os

import ROOT
import dd4hep
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Rectangle

ROOT.gROOT.SetBatch(True)

# volume-name prefixes defined by src/BarrelCalorimeterImaging_geo.cpp
LAYER = 'layer'
STAVE = 'stave'            # staves staveN<position> (-z) and staveP<position> (+z, rotated by 180 deg), air; readout id 'stave'
MODULE = 'AstroPix_Module'  # module envelope (air), readout id 'module', counted from the outer stave end
CHIP = 'AstroPix_Chip'     # chip assembly, readout id 'chip', counted from the outer end; daughter is the Si box

COLORS = ('#9ecae1', '#fdd0a2')  # even, odd module index
EDGE, MARK = '#333333', '#b2182b'


def find_volume(vol, name, depth=0):
    for i in range(vol.GetNdaughters()):
        dv = vol.GetNode(i).GetVolume()
        if dv.GetName() == name:
            return dv
        if depth < 4:
            r = find_volume(dv, name, depth + 1)
            if r:
                return r
    return None


def volid(node, field):
    '''readout id (physVolID) of a placement, e.g. module or chip'''
    return {k.decode() if isinstance(k, bytes) else k: v for k, v in dd4hep.PlacedVolume(node).volIDs()}[field]


def ty(node):
    '''y translation of a node in mm (only used to find the stave positions and order their staves)'''
    return 10 * node.GetMatrix().GetTranslation()[1]


def footprint(chain, shape):
    '''layer-local (y centre, y half-length, x half-width) in mm of a box shape placed through the node chain
    (outermost node first), using the full node matrices (staveP is rotated by 180 degrees)'''
    ys, xs = [], []
    for sx in (-1, 1):
        for sy in (-1, 1):
            p = array.array('d', [sx * shape.GetDX(), sy * shape.GetDY(), 0.])
            for n in reversed(chain):
                q = array.array('d', [0., 0., 0.])
                n.GetMatrix().LocalToMaster(p, q)
                p = q
            xs.append(10 * p[0])
            ys.append(10 * p[1])
    return dict(y=(max(ys) + min(ys)) / 2, hy=(max(ys) - min(ys)) / 2, hx=(max(xs) - min(xs)) / 2)


def extract(args):
    det = dd4hep.Detector.getInstance()
    det.fromXML(args.compact)
    sdet = det.detector(args.detector)
    if not sdet.isValid():
        raise SystemExit('detector %s not found in %s' % (args.detector, args.compact))
    sec = find_volume(sdet.placement().volume().ptr(), args.sector)
    if not sec:
        raise SystemExit('volume %s not found in %s' % (args.sector, args.detector))
    layer = [sec.GetNode(i) for i in range(sec.GetNdaughters())
             if sec.GetNode(i).GetName().startswith('%s%d_' % (LAYER, args.layer))]
    if not layer:
        raise SystemExit('layer %d not found in %s' % (args.layer, args.sector))
    lv = layer[0].GetVolume()
    # stave positions: staves sharing an x-z footprint, in placement order
    positions = {}
    for i in range(lv.GetNdaughters()):
        n = lv.GetNode(i)
        if n.GetName().startswith(STAVE):
            t = n.GetMatrix().GetTranslation()
            positions.setdefault((round(t[0], 4), round(t[2], 4)), []).append(n)
    positions = list(positions.values())
    if not 0 <= args.position < len(positions) or len(positions[args.position]) != 2:
        raise SystemExit('stave position %d: not found or not made of two staves' % args.position)
    staves = positions[args.position]
    staves.sort(key=ty)
    y0 = sum(ty(h) for h in staves) / 2  # stave position centre
    out = []
    for h in staves:
        hv = h.GetVolume()
        d = dict(id=volid(h, 'stave'), name=h.GetName(), mods=[], **footprint([h], hv.GetShape()))
        d['y'] -= y0
        for i in range(hv.GetNdaughters()):
            m = hv.GetNode(i)
            if not m.GetName().startswith(MODULE):
                continue
            mv = m.GetVolume()
            mod = dict(id=volid(m, 'module'), chips=[], **footprint([h, m], mv.GetShape()))
            mod['y'] -= y0
            for j in range(mv.GetNdaughters()):
                cn = mv.GetNode(j)
                c = cn.GetVolume().GetNode(0)  # Si box inside the chip assembly
                chip = dict(id=volid(cn, 'chip'), **footprint([h, m, cn, c], c.GetVolume().GetShape()))
                chip['y'] -= y0
                mod['chips'].append(chip)
            mod['chips'].sort(key=lambda c: c['y'])
            d['mods'].append(mod)
        d['mods'].sort(key=lambda m: m['y'])  # ascending y; the ids follow from the geometry
        out.append(d)
    return out


def box(ax, y, hy, hx, **kw):
    ax.add_patch(Rectangle((y - hy, -hx), 2 * hy, 2 * hx, **kw))


def draw(ax, staves, labels=False, lw=0.6):
    for h in staves:
        box(ax, h['y'], h['hy'], h['hx'], fill=False, ec=EDGE, ls='--', lw=1.4, zorder=4)
        for m in h['mods']:
            for c in m['chips']:
                box(ax, c['y'], c['hy'], c['hx'], fc=COLORS[m['id'] % 2], ec=EDGE, lw=lw, zorder=2)
            box(ax, m['y'], m['hy'], m['hx'], fill=False, ec=EDGE, lw=lw * 2, zorder=3)
            if labels:
                ax.text(m['y'], 0, str(m['id']), ha='center', va='center', fontsize=8, zorder=5)


def dim(ax, y1, y2, x, text, tx=None, ha='center', off=0.5):
    '''dimension marker for the gap y1..y2: two arrows pointing at its edges from outside, plus the text'''
    for a, b in ((y1 - off, y1), (y2 + off, y2)):
        ax.annotate('', xy=(b, x), xytext=(a, x), arrowprops=dict(arrowstyle='->', color=MARK, lw=1.3, shrinkA=0, shrinkB=0))
    ax.text((y1 + y2) / 2 if tx is None else tx, x + 0.35, text, color=MARK, ha=ha, va='bottom', fontsize=10, zorder=6,
            bbox=dict(fc='white', ec='none', alpha=0.85, pad=1))


def table(staves):
    h0, h1 = staves
    chip, mod = h0['mods'][0]['chips'][0], h0['mods'][0]
    cg = [b['y'] - b['hy'] - (a['y'] + a['hy']) for h in staves for m in h['mods'] for a, b in zip(m['chips'], m['chips'][1:])]
    mg = [b['y'] - b['hy'] - (a['y'] + a['hy']) for h in staves for a, b in zip(h['mods'], h['mods'][1:])]
    jg = h1['mods'][0]['y'] - h1['mods'][0]['hy'] - (h0['mods'][-1]['y'] + h0['mods'][-1]['hy'])
    eg = [(h['mods'][0]['y'] - h['mods'][0]['hy']) - (h['y'] - h['hy']) for h in staves[:1]] + \
         [(h['y'] + h['hy']) - (h['mods'][-1]['y'] + h['mods'][-1]['hy']) for h in staves[1:]]
    hg = (h1['y'] - h1['hy']) - (h0['y'] + h0['hy'])
    rows = [('chip Si size y x x [mm]', '%.3f x %.3f' % (2 * chip['hy'], 2 * chip['hx'])),
            ('chips per module', str(len(mod['chips']))),
            ('module envelope y x x [mm]', '%.3f x %.3f' % (2 * mod['hy'], 2 * mod['hx'])),
            ('stave envelope y x x [mm]', '%.3f x %.3f' % (2 * h0['hy'], 2 * h0['hx'])),
            ('modules per stave (staveN + staveP)', '%d + %d' % (len(h0['mods']), len(h1['mods']))),
            ('chips per stave position', str(sum(len(m['chips']) for h in staves for m in h['mods']))),
            ('chip-to-chip gap [mm] (min/max)', '%.3f / %.3f' % (min(cg), max(cg))),
            ('module-to-module gap [mm] (min/max)', '%.3f / %.3f' % (min(mg), max(mg))),
            ('junction module gap at y=0 [mm]', '%.3f' % jg),
            ('gap at the outer stave ends [mm]', '%.3f / %.3f' % tuple(eg)),
            ('staveN-staveP box gap at y=0 [mm]', '%.3f' % hg),
            ('stave position length (outer stave edges) [mm]', '%.3f' % (h1['y'] + h1['hy'] - (h0['y'] - h0['hy'])))]
    for k, v in rows:
        print('%-46s %s' % (k, v))
    print('%-46s %s' % ('module ids from -z to +z', ' '.join(str(m['id']) for h in staves for m in h['mods'])))
    print('%-46s %s' % ('chip ids in the first module of each stave',
                        ' | '.join('stave %d: ' % h['id'] + ' '.join(str(c['id']) for c in h['mods'][0]['chips'])
                                   for h in staves)))


def fig_stave(staves, args):
    h0, h1 = staves
    L = h1['y'] + h1['hy']
    fig = plt.figure(figsize=(16, 8.6))
    top = fig.add_axes([0.06, 0.64, 0.92, 0.25])
    draw(top, staves, labels=True)
    top.set_xlim(-L * 1.01, L * 1.01)
    top.set_ylim(-14, 14)
    top.set_yticks([])
    top.set_xlabel('y along the beam, from the stave centre [mm] (x across the stave is stretched)')
    for h in staves:
        ids = [m['id'] for m in (h['mods'] if h['y'] < 0 else h['mods'][::-1])]  # outer end first
        top.text(h['y'], 12, '%s (stave id %d, %sz): modules %d-%d from the outer end' %
                 (h['name'], h['id'], '-' if h['y'] < 0 else '+', ids[0], ids[-1]), ha='center', va='center', fontsize=10)
    # zoom targets: (title, y centre, annotations)
    m0 = h0['mods']
    gaps = [(h0['y'] - h0['hy'], 'outer stave end'),
            (0.5 * (m0[3]['chips'][4]['y'] - m0[3]['chips'][4]['hy'] + m0[3]['chips'][3]['y'] + m0[3]['chips'][3]['hy']),
             'chip-to-chip gap'),
            (0.5 * (m0[7]['y'] + m0[7]['hy'] + m0[8]['y'] - m0[8]['hy']), 'module-to-module gap'),
            (0.0, 'staveN-staveP junction (calorimeter centre)')]
    for k, (yc, title) in enumerate(gaps):
        ax = fig.add_axes([0.06 + 0.2375 * k, 0.07, 0.195, 0.40])
        draw(ax, staves, lw=0.8)
        w = 2.1
        ax.set_xlim(yc - w, yc + w)
        ax.set_ylim(10.4 - 2 * w, 10.4)
        ax.set_aspect('equal')
        ax.set_xlabel('y [mm]')
        ax.set_ylabel('x [mm]')
        ax.set_title('(%s) %s' % ('abcd'[k], title), fontsize=11)
        top.plot([yc, yc], [-10.2, 10.2], color=MARK, lw=1.5, zorder=6)
        top.text(yc, -12.3, '(%s)' % 'abcd'[k], color=MARK, ha='center', va='center', fontsize=11)
        fig.add_artist(ConnectionPatch((yc, -10.2), (0.5, 1), 'data', 'axes fraction', axesA=top, axesB=ax,
                                       color=MARK, lw=0.8, ls=':'))
        if k == 0:
            e = h0['y'] - h0['hy']
            dim(ax, e, m0[0]['y'] - m0[0]['hy'], 8.0, '%.3f mm' % (m0[0]['y'] - m0[0]['hy'] - e), tx=e + 1.2)
            ax.text(e - 0.1, 7.0, 'stave\nenvelope edge', ha='right', fontsize=9)
        elif k == 2:
            a, b = m0[7], m0[8]
            dim(ax, a['y'] + a['hy'], b['y'] - b['hy'], 8.0, '%.3f mm' % (b['y'] - b['hy'] - a['y'] - a['hy']))
        elif k == 1:
            a, b = m0[3]['chips'][3], m0[3]['chips'][4]
            dim(ax, a['y'] + a['hy'], b['y'] - b['hy'], 8.0, '%.3f mm' % (b['y'] - b['hy'] - a['y'] - a['hy']))
        else:
            a, b = m0[-1], h1['mods'][0]
            dim(ax, a['y'] + a['hy'], 0, 8.0, '%.3f mm' % (-a['y'] - a['hy']), tx=-0.3, ha='right', off=0.3)
            dim(ax, 0, b['y'] - b['hy'], 8.0, '%.3f mm' % (b['y'] - b['hy']), tx=0.3, ha='left', off=0.3)
            ax.text(0, 7.0, 'module gap %.3f mm\nstave envelopes touch' % (b['y'] - b['hy'] - a['y'] - a['hy']),
                    ha='center', fontsize=9, zorder=6, bbox=dict(fc='white', ec='none', alpha=0.9))
    fig.suptitle('BIC %s layer %d, stave position %d: %d chips (%s) - boxes: AstroPix chip Si (component0) footprints, '
                 'solid outlines: module envelopes (air), dashed: stave envelopes (air)' %
                 (args.sector, args.layer, args.position, sum(len(m['chips']) for h in staves for m in h['mods']),
                  'colour = module id parity'), fontsize=10.5, y=0.97)
    fig.savefig(args.output + '.png', dpi=150)


def fig_module(staves, args):
    fig, axes = plt.subplots(2, 1, figsize=(14, 5.6))
    for ax, h in zip(axes, staves):
        m = next(m for m in h['mods'] if m['id'] == 0)
        out = -1 if h['y'] < 0 else 1  # direction of the outer stave end along y, from the geometry
        for c in m['chips']:
            box(ax, c['y'], c['hy'], c['hx'], fc=COLORS[0], ec=EDGE, lw=0.8)
            ax.text(c['y'], 0, 'chip %d' % c['id'], ha='center', va='center', fontsize=9)
        box(ax, m['y'], m['hy'], m['hx'], fill=False, ec=EDGE, lw=1.6)
        c0 = m['chips'][0]
        ax.annotate('', xy=(c0['y'] - c0['hy'], c0['hx'] + 3), xytext=(c0['y'] + c0['hy'], c0['hx'] + 3),
                    arrowprops=dict(arrowstyle='<->', color=MARK, lw=1.3, shrinkA=0, shrinkB=0))
        ax.text(c0['y'], c0['hx'] + 4.5, 'chip %.3f mm' % (2 * c0['hy']), color=MARK, ha='center', fontsize=10)
        yb = -c0['hx'] - 3
        for sgn, txt in ((out, 'outer stave end'), (-out, 'towards junction (y = 0)')):
            ax.annotate(txt, xy=(m['y'] + sgn * (m['hy'] + 3.5), yb), xytext=(m['y'] + sgn * 12, yb), color=MARK,
                        ha='right' if sgn < 0 else 'left', va='center', fontsize=10,
                        arrowprops=dict(arrowstyle='->', color=MARK, lw=1.3, shrinkA=0, shrinkB=0))
        ax.set_xlim(m['y'] - m['hy'] - 4, m['y'] + m['hy'] + 4)
        ax.set_ylim(-c0['hx'] - 6, c0['hx'] + 8)
        ax.set_aspect('equal')
        ax.set_ylabel('x [mm]')
        ax.set_title('module %d of %s (stave id %d, %sz side), true scale: %.3f mm long, %d chips, '
                     '%.3f mm chip-to-chip gap (boxes: chip Si footprints, outline: module envelope)' %
                     (m['id'], h['name'], h['id'], '-' if h['y'] < 0 else '+', 2 * m['hy'], len(m['chips']),
                      m['chips'][1]['y'] - m['chips'][1]['hy'] - m['chips'][0]['y'] - m['chips'][0]['hy']), fontsize=10)
    axes[-1].set_xlabel('y from the stave centre [mm]')
    fig.tight_layout()
    fig.savefig(args.output + '_module.png', dpi=150)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--compact', help='compact detector file',
                        default=f'{os.getenv("DETECTOR_PATH")}/{os.getenv("DETECTOR_CONFIG")}.xml')
    parser.add_argument('-o', '--output', default='bic_stave_modules', help='output file prefix')
    parser.add_argument('--detector', default='EcalBarrelImaging', help='name of the imaging calorimeter')
    parser.add_argument('--sector', default='sector26', help='sector volume (sector1 for layer 1, sector26 for 2-6)')
    parser.add_argument('--layer', type=int, default=2, help='layer number')
    parser.add_argument('--position', type=int, default=0, help='stave position in the layer (staveN<n> and staveP<n>)')
    args = parser.parse_args()
    staves = extract(args)
    table(staves)
    fig_stave(staves, args)
    fig_module(staves, args)
