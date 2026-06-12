#!/bin/env python

import numpy as np
import os, sys, stat
import shutil
import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '../../')

import subprocess
import chromflock_common

resolution = 1e6;
folder = './'

print(f"Approximately {resolution} basepairs per bead")
mflock = chromflock_common.find_mflock()
chr_sizes = chromflock_common.chr_sizes

chr_sizes = np.concatenate([chr_sizes, chr_sizes])
nchr = len(chr_sizes)

print(f"{nchr} chromosomes")

nbeads = np.int64(chr_sizes/resolution);

labels = np.zeros(round(np.sum(nbeads)))

bstart = 0
for chr in range(0, nchr):
    labels[bstart:bstart+nbeads[chr]] = chr+1
    bstart = bstart + nbeads[chr]

# Adjust labels to match the color map (which can be changed)
# this step is not needed
labels[labels>24] = labels[labels>24]+32-24
np.save('labels.npy', labels.astype('uint8'))

shutil.copyfile('../../../src/mflock.lua', folder + '/mflock.lua')

cmd = [mflock, '--backbone', '-L', 'labels.npy',
       '--dconf', 'mflock.lua',
       '--outFolder', './',
       '--live',
       '--cmm']

cmds = " ".join(cmd)
print(f"command: {cmds}")
input("Press Enter to run...")
subprocess.run(cmd)
