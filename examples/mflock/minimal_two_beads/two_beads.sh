set -e

chromflock=../../../bin/chromflock
mflock=../../../bin/mflock
cp ../../../src/mflock.lua .

# The number of beads is specified by the label list
# The first chromosome is always 1 (0 is unused)
# Here we specify that there are two beads of chromosome 1
$chromflock string2any labels.u8 uint8_t 1 1

# Say that the first bead is connected to the second
# note that 0-indexing is used for the beads
$chromflock string2any contact_pairs.u32 uint32_t 0 1

$mflock --contact-pairs contact_pairs.u32 \
        --lFile labels.u8 \
        --outFolder './' \
        --config mflock.lua \
        --live
