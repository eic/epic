'''
    A script to visualize the fibers and of one grid from BEMC ScFi part
    Chao Peng (ANL)
'''
import os
import ROOT
import dd4hep
import DDRec
import argparse
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.ticker as ticker
from collections import OrderedDict


# helper function to do some type conversion for python wrapper of C++ function
def dict_to_cpp_vec(my_dict, dtype='int'):
    # use ROOT to make sure the type is correct
    vol_ids = ROOT.std.vector('pair<string, {}>'.format(dtype))()
    for field, fid in my_dict.items():
        vol_ids.push_back((field, fid))
    return vol_ids


# helper function to collect fibers under a grid
def get_fibers(detelem, id_desc, id_conv, id_dict, grid=1):
    gnode = detelem.volume().GetNode(grid-1)
    global_trans = detelem.nominal().worldTransformation()
    print(gnode.GetName())
    # get grid
    id_dict.update({'grid': grid, 'fiber': 0})
    gid = id_desc.encode(dict_to_cpp_vec(id_dict))
    # print(id_desc.decoder().valueString(vid))
    vol = id_conv.findContext(gid).volumePlacement().volume()
    grpos = id_conv.position(gid)
    # print(vol.boundingBox().dimensions())

    # use TGeoNode to get center positions
    # here it can also use converter to get position but it's much slower (adds an additional lookup process)
    fibers = []
    for i in np.arange(vol.GetNdaughters()):
        fnode = vol.GetNode(int(i))
        # NOTE, this is defined in geometry plugin, fiber_core is the wanted sensitive detector
        if 'fiber_core' not in fnode.GetName():
            continue
        fpos = np.array([0., 0., 0.])
        gpos = np.array([0., 0., 0.])
        pos = np.array([0., 0., 0.])
        fnode.LocalToMaster(np.array([0., 0., 0.]), fpos)
        gnode.LocalToMaster(fpos, gpos)
        # the two method below are equivalent
        global_trans.LocalToMaster(gpos, pos)
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
        fibers.append(pos)

    return fibers, grpos


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '-c', '--compact',
            dest='compact', required=True,
            help='Top-level xml file of the detector description.'
            )
    parser.add_argument(
            '--dpath', default='EcalBarrelScFi/stave0/layer1/slice1',
            dest='dpath',
            help='Path to the closest DetElement (slice) that contains the grids and fibers.'
            )
    parser.add_argument(
            '--outdir',
            dest='outdir', default='.',
            help='Output directory.'
            )
    parser.add_argument(
            '--readout',
            default='EcalBarrelScFiHits',
            help='Readout class for the EcalBarrelScFi.'
            )
    parser.add_argument(
            '--center-grid', dest='center_grid', type=int, default=2,
            help='The grid number for centering in the plot.'
            )
    args = parser.parse_args()

    # initialize dd4hep detector
    desc = dd4hep.Detector.getInstance()
    desc.fromXML(args.compact)
    converter = DDRec.CellIDPositionConverter(desc)

    # search the DetElement start from the world
    det_names = args.dpath.split('/')
    det = desc.world()

    for i, n in enumerate(det_names):
        children = det.children()
        try:
            det = det.child(n)
        except Exception:
            print('Failed to find DetElement {} from {}'.format(n, '/'.join(['world'] + det_names[:i])))
            print('Available children are listed below:')
            for n, d in children:
                print(' --- child: {}'.format(n))
            exit(-1)

    # build a volume manager so it triggers populating the volume IDs
    vman = dd4hep.VolumeManager(det, desc.readout(args.readout))

    # manager instances
    readout = desc.readout(args.readout)
    idspec = readout.idSpec()

    # build a dictionary for the readout ids
    mvol_id = det.volumeID()
    id_dict = OrderedDict()
    for idstr in idspec.decoder().valueString(mvol_id).split(','):
        field, fid = idstr.split(':')
        id_dict[field] = int(fid)

    # get radius of the fiber
    id_dict.update({'grid': args.center_grid, 'fiber': 1})
    fid = idspec.encode(dict_to_cpp_vec(id_dict))
    # need to do unit conversion from mm to cm (might from ROOT/DD4Hep difference)
    fr = converter.cellDimensions(fid)[0]/2./10.

    # plot fibers in the grid
    fig, ax = plt.subplots(figsize=(12, 12), dpi=160)
    # default color cycles
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # plot the center grid and adjacent grid
    grids = [
        (g, colors[i]) for i, g in enumerate([args.center_grid, args.center_grid - 1, args.center_grid + 1]) if g > 0
        ]
    for g, c in grids:
        # get position of fibers
        fibers, grid_pos = get_fibers(det, idspec, converter, id_dict, g)
        patches = []
        for fpos in fibers:
            patches.append(Circle((fpos[0], fpos[1]), fr))
        p = PatchCollection(patches, alpha=0.4, facecolors=(c,), edgecolors=('k',))
        ax.plot([grid_pos.X()], [grid_pos.Y()], marker='P', mfc=c, mec='k', ms=9, label='grid {} center'.format(g))
        ax.add_collection(p)
        if g == args.center_grid:
            ax.set_xlim(grid_pos.X()- 3., grid_pos.X() + 3.)
            ax.set_ylim(grid_pos.Y() - 3., grid_pos.Y() + 3.)
    ax.legend(fontsize=22)
    ax.tick_params(labelsize=20, direction='in')
    ax.set_xlabel('Global X (mm)', fontsize=22)
    ax.set_ylabel('Global Y (mm)', fontsize=22)
    ax.set_title('{} around grid {}'.format(args.dpath, args.center_grid), fontsize=22)

    fig.savefig(os.path.join(args.outdir, 'grid_fibers.png'))
