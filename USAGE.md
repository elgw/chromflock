# Usage:

If you have installed chromflock there will be man pages for the
binaries, as well a some short usage instructions via the **--help**
flag. The man pages, rendered to text, can also be viewed here:

 - [mflock](man/mflock.txt)
 - [aflock](man/aflock.txt)
 - [chromflock](man/chromflock.txt)
 - [cc2cpm](man/cc2cpm.txt)

## Minimal example
Shows how to generate 128 structures based on Hi-C data.

```
mkdir 128
cd 128
# Convert Hi-C/TCC data to a contact probability matrix and remove chrY if included
# Assuming that you already have some Hi-C data in binary format.
cc2cpm --hFile ../H.double --lFile ../L.uint8 --nStruct 8
chromflock init
# Edit settings
vim chromflock_gen
# Possibly also the dynamics settings
vim mflock.lua
# Then run
./chromflock_run
```

If `chromflock_run` was aborted for some reason, inspect `status.txt` to see the last thing that was finished and continue from any line, L, by
```
bash < (sed -n 'L,$p' chromflock_run)
```

The columns of the files `coords.csv` are `x`, `y`, `z`, `radius`, and
(if supplied) `preferred_radius`.


An example of 100 kb beads for a Haploid cell-line can be found in
[HAP1_100k.md](HAP1_100k.md).
