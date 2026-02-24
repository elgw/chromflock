# Shows how to use positional constraints or absolute positioning of
# beads via the --absolute command line argument.
#

import numpy as np
import subprocess

labels = np.array([1, 1, 1, 1, 1],
                  dtype=np.uint8)

# bead_id, x, y, z
absolute = np.array([[0, 1, 0, 0],
                     [4, 0, 0, 0]],
                    dtype=np.float32)
# no contacts
contact_pairs = np.array([],
                         dtype=np.uint32)

np.save('labels.npy', labels)
np.save('absolute.npy', absolute)
np.save('contact_pairs.npy', contact_pairs.astype(np.uint32))

cmd = ['../../../bin/mflock',
       '--absolute', 'absolute.npy',
       '--labels',  'labels.npy',
       '--config', '../../../src/mflock.lua',
       '--contact-pairs', 'contact_pairs.npy',
       '--outFolder', 'results',
       '--cmm']

ret = subprocess.run(cmd)

X = np.load('results/coords.npy')
for ab in absolute:
    pos = int(ab[0])
    got = X[pos, :-1]
    wanted = ab[1:]
    print(f'{got} {wanted}')
    assert(np.linalg.norm(got-wanted) < 1e-4)
