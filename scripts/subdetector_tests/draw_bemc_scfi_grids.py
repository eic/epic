'''
    A script to visualize the fibers and of one grid from BEMC ScFi part
    Chao Peng (ANL)
'''
import os
import ROOT
import dd4hep
import DDRec
import argparse
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from collections import OrderedDict


def dict_to_cpp_vec(my_dict, dtype='int'):
    # use ROOT to make sure the type is correct
    vol_ids = ROOT.std.vector('pair<string, {}>'.format(dtype))()
    for field, fid in my_dict.items():
        vol_ids.push_back((field, fid))
    return vol_ids



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
            help='Path to the scfi layer that contains the grids to draw.'
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
            '--center-grid', dest='center_grid', type=int, default=1,
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
    idspec = desc.readout(args.readout).idSpec()

    # build a dictionary for the readout ids
    mvol_id = det.volumeID()
    id_dict = OrderedDict()
    for idstr in idspec.decoder().valueString(mvol_id).split(','):
        field, fid = idstr.split(':')
        id_dict[field] = int(fid)

    # update
    id_dict.update({'grid': 1})
    vid = idspec.encode(dict_to_cpp_vec(id_dict))
    print(idspec.decoder().valueString(vid))
    vol = converter.findContext(vid).volumePlacement().volume()
    print(vol.GetNdaughters())


    # [(fid[0], int(fid[1])) for fid in idspec.decoder().valueString(mvol_id).split(',')]
    # use ROOT to make sure the type is correct
    # vol_ids = ROOT.std.vector('pair<string, int>')()
    #     vol_ids.push_back((field, int(fid)))
    # vol_ids[-2].second = 2
    # print(vol_ids)
    # print(idspec.encode(vol_ids))

    # dd4hep_decoder = description.readout(args.readout).idSpec().decoder()
    # lindex = dd4hep_decoder.index('x')
    # get_layer_id = np.vectorize(lambda cid: dd4hep_decoder.get(cid, lindex))
    # df.loc[:, 'layerID'] = get_layer_id(df['cellID'].astype(int).values)
