# Source code

```
aflock.c        # Assignment step
mflock.c        # Molecular dynamics step
  functional.c  # Defines the forces/vector field
  md.c          # Molecular dynamics schedule
  cmmwrite.c    # Write out Chimera cmm files
  liveview.c    # SDL visualization
    hsvrgb.c    # HSV -> RBG and RGB <- HSV color space conversion
wio.c           # Read/Write raw data using libz

string2any.c    # For testing, converts strings to raw data numbers
gen_cm.c        # For generating test data
```
