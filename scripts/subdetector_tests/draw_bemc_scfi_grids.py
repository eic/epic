# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Chao Peng
'''
    A script to visualize the fibers of some grids from BEMC ScFi part
    use case:
    python scripts/subdetector_tests/draw_bemc_scfi_grids.py -c epic_brycecanyon.xml
'''
import os
import ROOT
import dd4hep
import DDRec
import argparse
import numpy as np
from pydoc import locate
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.ticker as ticker
from collections import OrderedDict


# helper function to do some type conversion for python wrapper of C++ function
def dict_to_cpp_vec(my_dict, dtype='int'):
    # use ROOT to make sure the type is correct
    vol_ids = ROOT.std.vector('pair<string, {}>'.format(dtype))()
    dcast = locate(dtype)
    for field, fid in my_dict.items():
        vol_ids.push_back((field, dcast(fid)))
    return vol_ids


# helper function to collect fibers under a grid
def get_grid_fibers(det_elem, vol_man, id_conv, id_dict):
    # locate the nearest DetElement
    id_desc = vol_man.idSpec()
    try:
        # get fiber radius
        id_dict.update({'fiber': 1})
        fid = id_desc.encode(dict_to_cpp_vec(id_dict))
        # NOTE: for tube geometry, and it needs a cm to mm conversion
        fr = id_conv.cellDimensions(fid)[0]/2./10.

        # get the lowest level DetElement
        sdet = id_conv.findDetElement(id_conv.position(fid))
        gtrans = sdet.nominal().worldTransformation()

        # get grid node (it's not a DetElement)
        id_dict.update({'fiber': 0})
        gid = id_desc.encode(dict_to_cpp_vec(id_dict))
        gnode = id_conv.findContext(gid).volumePlacement()
        # print(id_desc.decoder().valueString(gid))
        grpos = id_conv.position(gid)
        grpos = np.array([grpos.X(), grpos.Y(), grpos.Z()])
    except Exception:
        return None, None

    # use TGeoNode to get center positions
    # here it can also use id_conv to do the same thing with cellIDs,
    # but it's much slower (adds an additional lookup process)
    fibers = []
    for i in np.arange(gnode.GetNdaughters()):
        fnode = gnode.GetDaughter(int(i))
        # NOTE, this is defined in geometry plugin, fiber_core is the only wanted sensitive detector
        if 'fiber_core' not in fnode.GetName():
            continue
        fpos = np.array([0., 0., 0.])
        gpos = np.array([0., 0., 0.])
        pos = np.array([0., 0., 0.])
        fnode.LocalToMaster(np.array([0., 0., 0.]), fpos)
        gnode.LocalToMaster(fpos, gpos)
        # the two method below are equivalent
        gtrans.LocalToMaster(gpos, pos)
        # detelem.nominal().localToWorld(gpos, pos)
        """ a test with converter
        if i < 50:
            id_dict.update({'fiber': int(len(fibers) + 1)})
            fid = id_desc.encode(dict_to_cpp_vec(id_dict))
            tpos = id_conv.position(fid)
            print(i,
                  fnode.GetName(),
                  np.asarray([tpos.X(), tpos.Y(), tpos.Z()]),
                  pos,
                  fpos,
                  gpos)
        """
        fibers.append(np.hstack([pos, fr]))

    return np.array(fibers), grpos


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '-c', '--compact',
            dest='compact', required=True,
            help='Top-level xml file of the detector description.'
            )
    parser.add_argument(
            '--detector', default='EcalBarrelScFi',
            dest='detector',
            help='Detector name.'
            )
    parser.add_argument(
            '--readout', default='EcalBarrelScFiHits',
            help='Readout class for the detector.'
            )
    parser.add_argument(
            '--grid-path', default='module:1,layer:6,slice:1,grid:3',
            help='Path down to a grid volume to be centered at the plot, with the format \"field1:i1,field2:i2,...\"',
            )
    parser.add_argument(
            '--outdir',
            dest='outdir', default='.',
            help='Output directory.'
            )
    parser.add_argument(
            '--adj-nlayers',
            dest='nlayers', type=int, default=2,
            help='number of adjacent layers to draw (+-n).'
            )
    parser.add_argument(
            '--adj-ngrids',
            dest='ngrids', type=int, default=2,
            help='number of adjacent grids to draw (+-n).'
            )
    parser.add_argument(
            '--window-size',
            dest='wsize', type=float, default=4.,
            help='Plot window size (mm).'
            )
    parser.add_argument(
            '--no-marker', action='store_true',
            help='Switch to draw a marker for grid center or not'
            )
    args = parser.parse_args()

    # initialize dd4hep detector
    desc = dd4hep.Detector.getInstance()
    desc.fromXML(args.compact)

    # search the detector
    det = desc.world()
    try:
        det = det.child(args.detector)
    except Exception:
        print('Failed to find detector {} from \"{}\"'.format(args.detector, args.compact))
        print('Available detectors are listed below:')
        for n, d in desc.world().children:
            print(' --- detector: {}'.format(n))
        exit(-1)

    # build a volume manager so it triggers populating the volume IDs
    vman = dd4hep.VolumeManager(det, desc.readout(args.readout))
    converter = DDRec.CellIDPositionConverter(desc)

    fields = OrderedDict([['system', det.id()]] + [v.split(':') for v in args.grid_path.split(',')])
    layer = int(fields.get('layer'))
    grid = int(fields.get('grid'))
    # add adjacent layers and grids, always put the central one (0, 0) first
    id_dicts = []
    for dl in np.hstack([np.arange(0, args.nlayers + 1), np.arange(-1, -args.nlayers - 1, step=-1)]):
        for dg in np.hstack([np.arange(0, args.ngrids + 1), np.arange(-1, -args.ngrids - 1, step=-1)]):
            if layer + dl < 1 or grid + dg < 1:
                continue
            new_dict = fields.copy()
            new_dict.update({'layer': layer + dl, 'grid': grid + dg})
            id_dicts.append(new_dict)

    # plot fibers in the grid
    fig, ax = plt.subplots(figsize=(12, 12), dpi=160)
    # default color cycle
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i, ids in enumerate(id_dicts):
        # color index number
        ic = (ids.get('grid') + (ids.get('layer') % 2)*4 - 1) % len(colors)
        c = colors[ic]
        fibers, gr_pos = get_grid_fibers(det, vman, converter, ids)
        if fibers is None:
            print('ignored {} because the volume might not exist.'.format(ids))
            continue

        patches = []
        for fi in fibers:
            patches.append(Circle((fi[0], fi[1]), fi[3]))
        p = PatchCollection(patches, alpha=0.6, facecolors=(c,), edgecolors=('k',))
        if not args.no_marker:
            ax.plot(gr_pos[0], gr_pos[1], marker='P', mfc=c, mec='k', ms=9, label='grid {}'.format(ids['grid']))
        ax.add_collection(p)
        # center at the first entry
        if i == 0:
            ax.set_xlim(gr_pos[0] - args.wsize, gr_pos[0] + args.wsize)
            ax.set_ylim(gr_pos[1] - args.wsize, gr_pos[1] + args.wsize)

    # ax.legend(fontsize=22)
    ax.tick_params(labelsize=20, direction='in')
    ax.set_xlabel('X (mm)', fontsize=22)
    ax.set_ylabel('Y (mm)', fontsize=22)
    ax.set_title('Centered at {}'.format('/'.join(['{}{}'.format(k, v) for k, v in fields.items()])), fontsize=22)
    fig.savefig(os.path.join(args.outdir, 'grid_fibers.png'))
