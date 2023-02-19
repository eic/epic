'''
    A script to visualize the fibers and of one grid from BEMC ScFi part
    Chao Peng (ANL)
'''
import os
import DDG4
import argparse
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '-c', '--compact',
            dest='compact', required=True,
            help='Top-level xml file of the detector description'
            )
    parser.add_argument(
            '-o',
            dest='outdir', default='.',
            help='Output directory.'
            )
    parser.add_argument(
            '--layer', type=int, default=1,
            help='The BEMC ScFi layer number.'
            )
    parser.add_argument(
            '--center-grid', dest='center_grid', type=int, default=1,
            help='The grid number for centering in the plot.'
            )
    args = parser.parse_args()

    # initialize dd4hep detector
    kernel = DDG4.Kernel()
    description = kernel.detectorDescription()
    kernel.loadGeometry("file:{}".format(args.compact))

    # dd4hep_decoder = description.readout(args.readout).idSpec().decoder()
    # lindex = dd4hep_decoder.index('x')
    # get_layer_id = np.vectorize(lambda cid: dd4hep_decoder.get(cid, lindex))
    # df.loc[:, 'layerID'] = get_layer_id(df['cellID'].astype(int).values)

    # always terminate dd4hep kernel
    kernel.terminate()

