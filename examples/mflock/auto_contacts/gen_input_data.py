import numpy as np

n = 160

L = np.zeros(n, dtype=np.uint8)
L = L + 1
L[0] = 0
L[n-1] = 2
np.save('labels.npy', L)

C = np.array([[0, 1], [n-2, n-1], [n-1, 0]], dtype=np.uint32)
np.save('contacts.npy', C)

"""
# The cmm file is not showing the auto-contacts

python ./gen_input_data.py
../../../build/mflock --backbone -L labels.npy --dconf mflock.lua --outFolder ./ --contact-pairs contacts.npy --vq 0.01 --cmm --autocontacts --live
chimerax cmmdump.cmm

"""
