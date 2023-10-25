Demonstrates how aflock and mflock can be used together.

``` shell
# Generate input data
./3beads_3structures.sh
# Run aflock to hand out initial contacts
export PATH=../../../bin:$PATH
./run_me.sh
# Run mflock on each structure
source mflock_jobs
```

This will generate the folders

``` shell
cf_000001
cf_000002
cf_000003
```

each containing a `cmmdump.cmm` file which can be visualized.

Please note that this isn't the intended way to use chromflock.
