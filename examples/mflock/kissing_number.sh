set -e

chromflock=../../bin/chromflock
mflock=../../bin/mflock
folder=kissing_number/

mkdir ${folder}

# 13 Beads in two colors
$chromflock string2any ${folder}labels.u8 uint8_t 1 1 1 1 1 1 1 1 1 1 1 1 2

# Connect the fist 12 to the 13th (note: 0-indexing)
$chromflock string2any ${folder}contact_pairs.u32 uint32_t 0 12 1 12 2 12 3 12 4 12 5 12 6 12 7 12 8 12 9 12 10 12 11 12

$mflock --contact-pairs ${folder}contact_pairs.u32 --lFile ${folder}labels.u8 --outFolder ${folder} --dconf ../../src/mflock.lua

# Note that the yellow connector in the structure indicates that bead
# 12 and 13 are adjacent (which makes more sense real data sets)
