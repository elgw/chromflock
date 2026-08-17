# Shows how to use positional constraints or absolute positioning of
# beads via the --absolute command line argument.
#
import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '../../')

import numpy as np
import subprocess
import chromflock_common

mflock = chromflock_common.find_mflock()
print(f'mflock: {mflock}')

labels = np.array([1, 1, 1, 1, 1],
                  dtype=np.uint8)

# bead_id, x, y, z
absolute = np.array([[0, 0.5, 0, 0],
                     [4, 0, 0, 0]],
                    dtype=np.float32)
# no contacts
contact_pairs = np.array([],
                         dtype=np.uint32)

np.save('labels.npy', labels)
np.save('absolute.npy', absolute)
np.save('contact_pairs.npy', contact_pairs.astype(np.uint32))

cmd = [mflock,
       '--absolute', 'absolute.npy',
       '--labels',  'labels.npy',
       '--config', '../../../src/mflock.lua',
       '--contact-pairs', 'contact_pairs.npy',
       '--outFolder', 'results',
       '--cmm']
print(f"cmd {' '.join(cmd)}")
ret = subprocess.run(cmd)

X = np.load('results/coords.npy')
print(X)
for ab in absolute:
    pos = int(ab[0])
    got = X[pos, :-1]
    wanted = ab[1:]
    print(f'got: {got} wanted: {wanted}', end='')
    if np.linalg.norm(got-wanted) > 0.1:
        print(f' !!! ')
    else:
        print(f' ok')
