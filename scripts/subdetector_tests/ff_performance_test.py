# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Chao Peng
'''
    A python script to benchmark the performance of far forward detectors with high energy proton/nuclei
    It tests each far forward detector (see FF_MOTHER and FF_COMP) and reads the process time at the end of the simulation

    Author: Chao Peng (ANL)
    Date: 07/19/2023
'''
import os
import re
import sys
import subprocess
import argparse
import shutil
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors


# script directory
SDIR = os.path.dirname(os.path.realpath(__file__))
# event generation script (particle gun is not working in this method, subprocess.run won't go through)
GEN_SCRIPT = os.path.join(SDIR, 'gen_particles.py')


def mrad2deg(th):
    return th/1000./np.pi*180.


def deg2mrad(ang):
    return ang/180.*np.pi*1000.


# the ff xml file that is used by the top-level xml file
# this script will modify the file so try to use a copy of the original file
FF_MOTHER = 'compact/far_forward_test.xml'
# ff components
# this list will be tested one-by-one (by modifying the mother xml file
FF_COMP = np.array([
        ('Beamline (Ion)', 'far_forward/ion_beamline.xml'),
        ('Beamline ($e^-$)', 'far_forward/electron_beamline.xml'),
        ('B0 Beampipe', 'far_forward/beampipe_hadron_B0.xml'),
        ('B0 Tracker', 'far_forward/B0_tracker.xml'),
        ('B0 ECal', 'far_forward/B0_ECal.xml'),
        ('Off-M Tracker', 'offM_tracker.xml'),
        ('ZDC', 'far_forward/ZDC.xml'),
        ('Roman Pots', 'far_forward/roman_pots_eRD24_design.xml'),
])


# kwargs are used by gen_cmd and sim_cmd
def sim_performance_test(**kwargs):
    gen_file = 'ff_test_gen.hepmc'
    # generate particles
    gen_cmd = [
        'python', GEN_SCRIPT, gen_file,
        '-n={nev}',
        '--angmin={angle_min}',
        '--angmax={angle_max}',
        '--pmin={p_min}',
        '--pmax={p_max}',
        '--phmin={phi_min}',
        '--phmax={phi_max}',
        '--particles={particles}',
        ]
    gen_cmd = [c.format(**kwargs) for c in gen_cmd]
    print(' '.join(gen_cmd))
    subprocess.run(gen_cmd, check=True)

    # simulation
    sim_cmd = [
        'ddsim',
        '--runType=batch',
        '--part.minimalKineticEnergy=1*TeV',
        '--filter.tracker=edep0',
        # '-v=WARNING',
        '--numberOfEvents={nev}',
        # '--physics.list {physics_list}',
        # '--enableGun',
        # '--gun.particle=proton', '--gun.energy=275*GeV',
        # '--gun.position=\"0.0 0.0 0.0*cm\"',
        # '--gun.phiMin=0.', '--gun.phiMax=2.*pi',
        # '--gun.thetaMin=0.003', '--gun.thetaMax=0.004',
        # '--gun.distribution=uniform',
        '--inputFiles={}'.format(gen_file),
        '--outputFile={sim_file}',
        '--compact={compact}',
        ]
    sim_cmd = [c.format(**kwargs) for c in sim_cmd]
    print(' '.join(sim_cmd))
    p = subprocess.run(sim_cmd, stdout=subprocess.PIPE)
    lines = p.stdout.decode('utf-8').split('\n')[-4:]

    pat = re.compile('Event Processing:\s+([-+]?(?:\d*\.*\d+))\s+s')
    r = pat.search('\n'.join(lines), re.MULTILINE)
    if not r:
        print('Cannot find event processing time in:')
        print('\n'.join(lines))
        return None
    return float(r.group(1))


if __name__ == '__main__':
    # argument parser
    parser = argparse.ArgumentParser()

    parser.add_argument('compact',
            help='A Top-level XML file of the detector discription.'
            )
    parser.add_argument('--ff-compact',
            default=FF_MOTHER,
            help='The XML file for far-forward detectors, used by the top-level XML file. It will be modified by the script.'
            )
    parser.add_argument('--nev',
            default=100, type=int,
            help='Number of events.'
            )
    parser.add_argument('--particles',
            default='proton',
            help='Type of particles.'
            )
    parser.add_argument('--sim-file',
            default='ff_test.edm4hep.root',
            help='Temporary output root file.'
            )
    parser.add_argument('--theta-min',
            default=3, type=float,
            help='Min. polar angle (mrad).'
            )
    parser.add_argument('--theta-max',
            default=3.5, type=float,
            help='Max. polar angle (mrad).'
            )
    parser.add_argument('--p-min',
            default=275, type=float,
            help='Min. momentum of the particles (GeV).'
            )
    parser.add_argument('--p-max',
            default=275, type=float,
            help='Max. momentum of the particles (GeV).'
            )
    parser.add_argument('--phi-min',
            default=0, type=float,
            help='Min. phi angle of the particles (degree).'
            )
    parser.add_argument('--phi-max',
            default=360, type=float,
            help='Max. phi angle of the particles (degree).'
            )

    args = parser.parse_args()
    kwargs = vars(args)
    # convert mrad to angle
    kwargs['angle_min'] = mrad2deg(args.theta_min)
    kwargs['angle_max'] = mrad2deg(args.theta_max)

    # all ff components
    t = sim_performance_test(**kwargs)
    result = [('All FF', t)]

    # read ff-compact
    file1 = open(args.ff_compact, 'r')
    lines = file1.readlines()
    file1.close()

    for comp1, xml1 in FF_COMP:
        print('Running test for ff component: {}'.format(comp1))
        new_lines = []
        for nl in lines:
            is_needed = True
            for comp2, xml2 in FF_COMP:
                if xml2 in nl and comp2 != comp1:
                    is_needed = False
            if is_needed:
                new_lines.append(nl)
        file2 = open(args.ff_compact, 'w')
        file2.writelines(new_lines)
        file2.close()

        # print(comp1)
        # print(''.join(lines))
        # print(''.join(new_lines))

        t = sim_performance_test(**vars(args))
        result.append((comp1, t))
        print('{}: {:.2f} s for {:d} events'.format(comp1, t, 100))

    # restore it back
    file2 = open(args.ff_compact, 'w')
    file2.writelines(lines)
    file2.close()


    # Build the plot
    result = np.array(result)
    p = (args.p_min + args.p_max)/2.    # GeV
    theta = (args.theta_min + args.theta_max)/2.    # mrad
    nev = float(args.nev)

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    fig, ax = plt.subplots(figsize=(8, 6), dpi=160,  gridspec_kw=dict(top=0.95, bottom=0.25, left=0.15, right=0.98))
    x_pos = np.arange(len(result))
    for i, (x, d, e) in enumerate(zip(x_pos, result.T[1].astype(float)/nev, result.T[1].astype(float)/np.sqrt(nev)/nev)):
        ic = colors[i % len(colors)]
        ax.bar(x, d, yerr=e, align='center', color=mcolors.to_rgba(ic, alpha=0.5), ec=ic, ecolor='black', capsize=10)
    ax.set_ylabel('Process Time (s / event)', fontsize=16)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(result.T[0], rotation=45, ha='right', fontsize=14)
    ax.set_title('{:.0f} GeV/c {} @ ~{:.1f} mrad'.format(p, args.particles, theta), fontsize=16)
    ax.yaxis.grid(ls=':')
    ax.tick_params(labelsize=16)
    ax.set_axisbelow(True)

    # Save the figure
    # fig.tight_layout()
    fig.savefig('ff_test_{:.0f}_{}_{:.0f}_mrad.png'.format(p, args.particles, theta))
