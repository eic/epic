# produce auxiliary ROOT file for the IRT algorithm

import shutil, os, sys, argparse
import xml.etree.ElementTree as et

parser = argparse.ArgumentParser()
parser.add_argument(
        '-o', '--output-name', dest='outFile', default='irt-drich.root',
        help='output auxiliary ROOT file name', type=str
        )
argv = parser.parse_args()

import DDG4


def run():

    # compact files
    if not 'DETECTOR_PATH' in os.environ:
        print('ERROR: env var DETECTOR_PATH not set',file=sys.stderr)
        exit(1)
    mainFile = os.environ['DETECTOR_PATH'] + '/ecce.xml'
    richFile = os.environ['DETECTOR_PATH'] + '/compact/drich.xml'

    # backup original richFile, then parse
    shutil.copy(richFile,richFile+'.bak')
    richTree = et.parse(richFile)

    # enable `DRICH_create_irt_file` mode
    for constant in richTree.iter(tag='constant'):
        if(constant.attrib['name']=='DRICH_create_irt_file'):
            constant.set('value','1')

    # set auxiliary file name
    for detector in richTree.iter(tag='detector'):
        detector.set('irt_filename',argv.outFile)

    # overwrite original richFile
    richTree.write(richFile)

    # produce IRT config file
    try:
        kernel = DDG4.Kernel()
        kernel.loadGeometry(f'file:{mainFile}')
        kernel.terminate()
        print(f'\n -> produced {argv.outFile}\n')
    except:
        pass

    # revert to the original richFile
    os.replace(richFile+'.bak',richFile)

if __name__ == "__main__":
    run()
