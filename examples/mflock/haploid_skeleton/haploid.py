#!/bin/env python

import numpy as np
import os, sys, stat
import shutil

resolution = 1e6;
folder = f'./'
labels_out = 'labels.npy'
pairs_out = 'pairs.npy'

print(f"Approximately {resolution} basepairs per bead")


chr_sizes = np.array([247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754])
nchr = len(chr_sizes)


nbeads = np.int64(chr_sizes/resolution);

labels = np.zeros(round(np.sum(nbeads)))

bstart = 0
for chr in range(0, nchr):
    labels[bstart:bstart+nbeads[chr]] = chr+1
    bstart = bstart + nbeads[chr]


#labels.astype('uint8').tofile(folder + '/labels.u8')
np.save(labels_out, labels.astype('uint8'))


# Create the backbone contact pairs
P = np.zeros([labels.shape[0], 2], np.uint32)
idx = 0;
for bead in range(0, len(labels)-1):
    if(labels[bead] == labels[bead+1]):
        P[idx, :] = [bead, bead+1]
        idx = idx + 1
P = P[0:idx, :]


# P.astype('uint32').tofile(folder + '/contact_pairs.u32')
np.save(pairs_out, P)

shutil.copyfile('../../..//src/mflock.lua', folder + '/mflock.lua')

print("Run the follow command to continue:")
print("")
print(f"mflock --contact-pairs {pairs_out} -L {labels_out} --dconf mflock.lua --outFolder ./ --live --cmm\n");
