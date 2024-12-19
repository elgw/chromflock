set -e

chromflock=../../../bin/chromflock
mflock=../../../bin/mflock

if [ ! -f "mflock.lua" ]; then
    cp ../../../src/mflock.lua .
fi


# 10 beads
$chromflock string2any labels.u8 uint8_t \
            1 2 2 2 2 \
            2 2 2 2 3

# Link up the beads
$chromflock string2any contact_pairs.u32 uint32_t \
            0 1 \
            1 2 \
            2 3 \
            3 4 \
            4 5 \
            5 6 \
            6 7 \
            7 8 \
            8 9 \

# Tell bead i0 to be to the right
# Tell bead i9 to be to the left
$chromflock string2any bead_wells.f64 double \
            0 0  0.9 0 \
            0 0 -0.9 0 \
            9 0 0 0 \

$mflock --contact-pairs contact_pairs.u32 \
        --lFile labels.u8 \
        --outFolder ./ \
        --dconf mflock.lua \
        --bead-wells bead_wells.f64 \
        --verbose 2 \
        --vq 0.01 \
        --live
