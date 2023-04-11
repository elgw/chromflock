#!/bin/bash
set -e

NAME="Erik Wernersson"
EMAIL="https://github.com/elgw/chromflock/issues"

set -e # exit on error

if [ $EUID = 0 ]; then
    echo "WARNING: You are running as root, please abort!"
    echo "Sleeping for 10 s"
    sleep 10
fi

savedir="$(pwd)"
builddir=$(mktemp -d)
cd $builddir
pkgdir=$(mktemp -d)/chromflock
echo "pkgdir=$pkgdir"

ver_major=`sed -rn 's/^#define.*CF_VERSION_MAJOR.*([0-9]+).*$/\1/p' < $savedir/src/cf_version.h`
ver_minor=`sed -rn 's/^#define.*CF_VERSION_MINOR.*([0-9]+).*$/\1/p' < $savedir/src/cf_version.h`
ver_patch=`sed -rn 's/^#define.*CF_VERSION_PATCH.*([0-9]+).*$/\1/p' < $savedir/src/cf_version.h`
pkgver="${ver_major}.${ver_minor}.${ver_patch}"
echo "pkgver=$pkgver"
arch=$(dpkg --print-architecture)
echo "arch=$arch"

echo "Preparing files"
# Copy files to the correct places
# binaries
mkdir -p $pkgdir/usr/local/bin
cp $savedir/bin/mflock $pkgdir/usr/local/bin/mflock
cp $savedir/bin/aflock $pkgdir/usr/local/bin/aflock
cp $savedir/bin/chromflock $pkgdir/usr/local/bin/chromflock

# man pages
mkdir -p $pkgdir/usr/share/man/man1/
cp $savedir/man/mflock.1 $pkgdir/usr/share/man/man1/
cp $savedir/man/aflock.1 $pkgdir/usr/share/man/man1/
cp $savedir/man/chromflock.1 $pkgdir/usr/share/man/man1/

# data files
datadir=$pkgdir/usr/local/share/chromflock/
mkdir -p $datadir
cp $savedir/src/mflock.lua $datadir/
cp $savedir/util/chromflock_gen $datadir/

cd $pkgdir
size=$(du -k usr | tail -n1 | sed 's/usr//')

mkdir $pkgdir/DEBIAN/
cat > $pkgdir/DEBIAN/control << END
Package: chromflock
Version: $pkgver
Architecture: $arch
Depends: libcairo2, libreadline8, zlib1g, libsdl2-2.0-0, parallel
Conflicts:
Maintainer: $NAME <$EMAIL>
Installed-Size: $size
Section: custom
Priority: optional
Homepage: https://github.com/elgw/chromflock/
Description: Chromflock
END

cd $pkgdir/../
dpkg-deb --build chromflock chromflock_${pkgver}_${arch}.deb
mv chromflock_${pkgver}_${arch}.deb "$savedir"
