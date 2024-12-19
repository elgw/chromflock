# CHANGELOG

## 0.4.3
- Fixed reading and writing empty contact lists.

## 0.4.2
- Fixed broken `--diploid` mode.

## 0.4.1
- Added `--bead-wells` to **mflock** to add bead-specific attractors.
- Added `--diploid` to **mflock** to explicitly inform that the labels
  should be duplicated (since diploid structures are generated from a
  Hi-C matrix which does not differentiate between the two copies).
- Added some toy examples under `/examples`
- Added the color map under `/doc`

Broken:
- `make SDL=0`.

## 0.4.0

- Sparse representation of contact restraints Lists of contact pairs,
  replace the contact indicator matrices previously used. This change
  does significantly reduce the storage required for a dataset and
  also reduces the processing times somewhat. On the other hand it
  breaks compatibility with older versions.

  The overall usage procedure should not be changed though. The typical
  chain of commands:

  ``` shell
  chromflock init
  edit chromflock_gen
  ./chromflock_gen
  ./chromflock_run
  ```
  is the same.

- Refactoring and commenting: Parts of the code has been refactored to
be easier to read and understand. Comments has been added along the
way but the source is still far from being fully documented.

makefile:
- The `-fanalyzer` flag has been added to the debug builds. Currently
  no warnings are issued for any of the binaries with GCC 11.4.

aflock:
- Removed unused command line arguments.
- Output files have new names and format. `afock -F` does now produce:
  `assigned_contacts.u16` and `measured_contacts.u16`.
- Assigned contacts are stored as `contact-pairs.u32.gz` in the
  structure folders.
- RAM requirements reduced to approximately 1/4 of the previous
  version due to change of data types.

mflock:
- Command line arguments changed.
- Removed the option to set the simulation parameters from command line.
- Reads/writes gz compressed contact pairs instead of the binary
  indication matrix (W.uint8).
- Writes ISO time to log files at start and finish.
- Writes the command line to the log file.

chromflock string2any/any2string
- Added `uint16_t` and `uint32_t` to the list of supported data
  types. Also short aliases will work, i.e.: `u8`, `u16`, `u32`.

## 0.3.9
- Compiles with **-pedantic** and **-fanalyzer** without any warnings
  with gcc 11.4.0.

## 0.3.8
- Fixed issues with the makefile. SDL is used by default.
- Automatic RAM usage limitation (only on Linux) which should give a
  more graceful exit than **SEGMENTATION FAULT** when running out of
  RAM. The process is quite crude so it might not work in all
  situations (only checks the amount of free RAM once).

## 0.3.7
- Fixed infinite loop when finding scaling factor in hic2cpm (again!).
- Fixed the **--aOut** argument which was ignored (hic2cpm).

## 0.3.6
- Added the `--mean` option to `chromflock hic2cpm`.

## 0.3.5
- Fixed issues that `hic2cpm` sometimes didn't find a scaling value.

## 0.3.4
- Cleaned up the makefile
- `mflock.lua` is now created by `chromflock init`, forgot that in the
  previous version.

## 0.3.3
- New behaviour, instead of `chromflock`, do `chromflock init`.
- `cc2cpm` is baked into the binary `chromflock` and renamed to
  `hic2cpm`
- `sprite2cpm`, `any2string` and `string2any` also moved to `chromflock`.
- documentation is based on markdown instead of troff.

## 0.3.2
 - Cleaned up the makefile.
 - Added **--version** argument to **aflock**, **mflock**, and **cc2cpm**
 - Fixed typos in the **cc2cpm** documentation.

## 0.3.1
 - Added script to make a deb package for Ubuntu. That is the
   preferred way to install from now.
