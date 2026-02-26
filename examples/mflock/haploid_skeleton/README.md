Test mflock with chromosome skeletons at 1 Mb resolution or 3010 beads.

Adjacent beads are connected to form chromosomes but there are no
further contacts.

``` shell
# Generate input data
./hapoloid
mflock --contact-pairs contact_pairs.npy -L labels.npy --dconf mflock.lua --outFolder ./ --live --cmm
```
