"""
Demonstrate how a custom color map is used

In this case, there will be 3 beads with label 50, 100, 150
which will be colored red, green and blue

Generate the cmm file with
$ ../../../bin/mflock --labels labels.npy --cmap cmap.npy --cmm  --config ../../../src/mflock.lua --outFolder ./

and then open with

$ chimerax cmmdump.cmm

"""

import numpy as np

cmap = np.zeros((256, 3))

cmap[50, :] = [255, 0, 0]
cmap[100, :] = [0, 255, 0]
cmap[150, :] = [0, 0, 255]

np.save("cmap.npy", cmap.astype('uint8'))

labels = np.array([50, 100, 150])
np.save('labels.npy', labels)
