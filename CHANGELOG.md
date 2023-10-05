# CHANGELOG

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
