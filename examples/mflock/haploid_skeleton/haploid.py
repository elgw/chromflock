#!/bin/env python

import numpy as np
import os, sys, stat
import shutil

import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '../../')

import numpy as np
import subprocess
import chromflock_common as cf

mflock = cf.find_mflock()
chr_sizes = cf.chr_sizes

resolution = 1e6;
folder = f'./'
labels_out = 'labels.npy'
pairs_out = 'pairs.npy'

print(f"Approximately {resolution} basepairs per bead")

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
print(f"{mflock} --contact-pairs {pairs_out} -L {labels_out} --dconf mflock.lua --outFolder ./ --live --cmm\n");
print("")
print("should give the same results as")
print(f"{mflock} --backbone -L {labels_out} --dconf mflock.lua --outFolder ./ --live --cmm\n");
