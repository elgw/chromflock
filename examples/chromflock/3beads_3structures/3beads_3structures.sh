set -e

nstruct=3
chromflock=../../../bin/chromflock
mflock=../../../bin/mflock
aflock=../../../bin/aflock

cp ../../../src/mflock.lua ./

# The number of beads is specified by the label list
# The first chromosome is always 1 (0 is unused)
# Here we specify that there are two beads of chromosome 1
$chromflock string2any labels.u8 uint8_t 1 1 1

# Say that the first bead is connected to the second
# note that 0-indexing is used for the beads
$chromflock string2any contact_probability_matrix.double double \
            0.0 1.0 0.4 \
            1.0 0.0 1.0 \
            0.4 1.0 0.0 \

echo "# This will initiate ${nstruct} structures"
echo "aflock --A contact_probability_matrix.double --nStruct ${nstruct} --init --vq 0.2 --mArgs \"--dconf mflock.lua --labels labels.u8\"" > run_me.sh
chmod +x run_me.sh

echo "Everything set up"
echo "Continue with:"
echo "\$ ./run_me.sh"
