import numpy as np

"""

First example,

label 0: bead #0
label 1: beads #1 -- #n-2
label 2: bead #n-1

Using the --backbone argument connects beads of the same label (linearly) so
the beads with label 1 will be connected to form a string.

Then we add contacts between bead #0 and #1, #n-2, #n-1 and $n-1 and 0. Using the
--autocontacts argument those contacts will only be enabled when the corresponding beads are in
proximity.

It is possible that all beads connect to form a circle. But in most
cases at least one of the single beads finds its way to the string of
label 1.

Please note that the --live does not show any of the contacts.
Also in the cmm file the enabled auto contacts are not shown (fix this!).

Run with:

python ./gen_input_data.py
../../../build/mflock --backbone -L labels.npy --dconf mflock.lua --outFolder ./ --contact-pairs contacts.npy --vq 0.01 --cmm --autocontacts --live
chimerax cmmdump.cmm


../../../build/mflock --backbone -L labels2.npy --dconf mflock.lua --outFolder ./ --contact-pairs contacts2.npy --vq 0.01 --cmm --autocontacts --live
chimerax cmmdump.cmm


"""
n = 160

L = np.zeros(n, dtype=np.uint8)
L = L + 1
L[0] = 0
L[n-1] = 2
np.save('labels.npy', L)

C = np.array([[0, 1], [0, n-1], [n-2, n-1]], dtype=np.uint32)
np.save('contacts.npy', C)


# Second example, autocontacts to connect a lot of things

L = 0*L
np.save('labels2.npy', L)
C = np.zeros( [int(n*(n-1)/2), 2], dtype=np.uint32)
pos = 0
for i in range(0, n):
    for j in range(i+1, n):
        C[pos, :] = [i, j]
        pos+=1
np.save('contacts2.npy', C)
