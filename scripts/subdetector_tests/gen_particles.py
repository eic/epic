# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Chao Peng
'''
    A simple script to generate single particles in HEPMC3 format

    Author: Chao Peng (ANL)
    Date: 07/19/2023
'''
import os
import sys
from pyHepMC3 import HepMC3 as hm
import numpy as np
import argparse


PARTICLES = {
    "pion0": (111, 0.1349766),       # pi0
    "pion+": (211, 0.13957018),      # pi+
    "pion-": (-211, 0.13957018),     # pi-
    "kaon0": (311, 0.497648),        # K0
    "kaon+": (321, 0.493677),        # K+
    "kaon-": (-321, 0.493677),       # K-
    "proton": (2212, 0.938272),      # proton
    "neutron": (2112, 0.939565),     # neutron
    "electron": (11, 0.51099895e-3), # electron
    "positron": (-11, 0.51099895e-3),# positron
    "photon": (22, 0),               # photon
    "muon": (13, 105.6583755),       # muon
}


def gen_event(p, theta, phi, pid, mass):
    evt = hm.GenEvent(hm.Units.MomentumUnit.GEV, hm.Units.LengthUnit.MM)
    # final state
    state = 1
    e0 = np.sqrt(p*p + mass*mass)
    px = np.cos(phi)*np.sin(theta)
    py = np.sin(phi)*np.sin(theta)
    pz = np.cos(theta)

    # beam
    pbeam = hm.GenParticle(hm.FourVector(0, 0, 0, 0.938272), 2212, 4)
    ebeam = hm.GenParticle(hm.FourVector(0, 0, e0, np.sqrt(e0*e0 + 0.511e-3*0.511e-3)), 11, 4)

    # out particle
    hout = hm.GenParticle(hm.FourVector(px*p, py*p, pz*p, e0), pid, state)

    # vertex
    vert = hm.GenVertex()
    vert.add_particle_in(ebeam)
    vert.add_particle_in(pbeam)
    vert.add_particle_out(hout)
    evt.add_vertex(vert)
    return evt


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('output', help='path to the output file')
    parser.add_argument('-n', type=int, default=1000, dest='nev', help='number of events to generate')
    parser.add_argument('-s', type=int, default=-1, dest='seed', help='seed for random generator')
    parser.add_argument('--parray', type=str, default="", dest='parray',
                        help='an array of momenta in GeV, separated by \",\"')
    parser.add_argument('--pmin', type=float, default=8.0, dest='pmin', help='minimum momentum in GeV')
    parser.add_argument('--pmax', type=float, default=100.0, dest='pmax', help='maximum momentum in GeV')
    parser.add_argument('--angmin', type=float, default=0.0, dest='angmin', help='minimum angle in degree')
    parser.add_argument('--angmax', type=float, default=20.0, dest='angmax', help='maximum angle in degree')
    parser.add_argument('--phmin', type=float, default=0.0, dest='phmin', help='minimum angle in degree')
    parser.add_argument('--phmax', type=float, default=360.0, dest='phmax', help='maximum angle in degree')
    parser.add_argument('--particles', type=str, default='electron', dest='particles',
                        help='particle names, support {}'.format(list(PARTICLES.keys())))

    args = parser.parse_args()

    # random seed (< 0 will get it from enviroment variable 'SEED', or a system random number)
    if args.seed < 0:
        args.seed = os.environ.get('SEED', int.from_bytes(os.urandom(4), byteorder='big', signed=False))
    print("Random seed is {}".format(args.seed))
    np.random.seed(args.seed)

    output = hm.WriterAscii(args.output);
    if output.failed():
        print("Cannot open file \"{}\"".format(args.output))
        sys.exit(2)

    # build particle info
    parts = []
    for pid in args.particles.split(','):
        pid = pid.strip()
        if pid not in PARTICLES.keys():
            print('pid {:d} not found in dictionary, ignored.'.format(pid))
            continue
        parts.append(PARTICLES[pid])

    # p values
    pvals = np.random.uniform(args.pmin, args.pmax, args.nev) if not args.parray else \
            np.random.choice([float(p.strip()) for p in args.parray.split(',')], args.nev)
    thvals = np.random.uniform(args.angmin, args.angmax, args.nev)/180.*np.pi
    phivals = np.random.uniform(args.phmin, args.phmax, args.nev)/180.*np.pi
    partvals = [parts[i] for i in np.random.choice(len(parts), args.nev)]

    count = 0
    for p, theta, phi, (pid, mass) in zip(pvals, thvals, phivals, partvals):
        if (count % 1000 == 0):
            print("Generated {} events".format(count), end='\r')
        evt = gen_event(p, theta, phi, pid, mass)
        output.write_event(evt)
        evt.clear()
        count += 1

    print("Generated {} events".format(args.nev))
    output.close()
