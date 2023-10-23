#!/bin/env python

import numpy as np
import os, sys, stat
import shutil

resolution = 1e6;
folder = f'demo_haploid_{resolution:.1e}'

print(f"Approximately {resolution} basepairs per bead")
print(f"Generating folder: {folder}")
os.mkdir(folder)

chr_sizes = np.array([247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754])
nchr = len(chr_sizes)


nbeads = np.int64(chr_sizes/resolution);

labels = np.zeros(round(np.sum(nbeads)))

bstart = 0
for chr in range(0, nchr):
    labels[bstart:bstart+nbeads[chr]] = chr+1
    bstart = bstart + nbeads[chr]


labels.astype('uint8').tofile(folder + '/labels.u8')

P = []
# Create the backbone contact pairs
for bead in range(0, len(labels)-1):
    if(labels[bead] == labels[bead+1]):
        P = np.append(P, [bead, bead+1])


P.astype('uint32').tofile(folder + '/contact_pairs.u32')

shutil.copyfile('../../src/mflock.lua', folder + '/mflock.lua')

with open(folder + '/run_me.sh', 'w') as fid:
    fid.write('set -e\n')
    fid.write('../../../bin/mflock --contact-pairs contact_pairs.u32 -L labels.u8 --dconf mflock.lua --outFolder ./\n');

os.chmod(folder + '/run_me.sh', stat.S_IRWXU)

print("Done")
