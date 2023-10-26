Test mflock with chromosome skeletons at 1 Mb resolution or 2x3010 beads.

Adjacent beads are connected to form chromosomes but there are no
further contacts.

In this case the `--diploid` flag is used to chromflock which will
duplicate the labels saved in `labels.u8` for the second copy.

``` shell
# Generate input data
./haploid
# Run mflock with --live
./run_me.sh
```
