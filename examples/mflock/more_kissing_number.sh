set -e

chromflock=../../bin/chromflock
mflock=../../bin/mflock
folder=more_kissing_number/

mkdir ${folder}

# 20 Beads in two colors
# 19 connections to the last bead
$chromflock string2any ${folder}labels.u8 uint8_t 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2

# Connect the fist 12 to the 13th (note: 0-indexing)
$chromflock string2any ${folder}contact_pairs.u32 uint32_t 0 19 1 19 2 19 3 19 4 19 5 19 6 19 7 19 8 19 9 19 10 19 11 19 12 19 13 19 14 19 15 19 16 19 17 19 18 19

$mflock --contact-pairs ${folder}contact_pairs.u32 --lFile ${folder}labels.u8 --outFolder ${folder} --dconf ../../src/mflock.lua

# Note that the yellow connector in the structure indicates that bead
# 12 and 13 are adjacent (which makes more sense real data sets)
