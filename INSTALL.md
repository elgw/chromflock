# INSTALLATION

On Ubuntu 22.04

``` shell
cd src/lua-5.3.5/
make linux # or pick another suitable target
cd ../../
make
./makedeb-ubuntu_2204
sudo apt install ./chromflock_*.deb
```

Check the dependencies section if there are any missing libraries. If
you have another OS, please check the section **Custom install**.

## Dependencies

For compilation the following is needed:

 * gcc
 * libreadline -- to compile lua
 * zlib -- to write `.gz` files
 * pkg-config -- for the makefile
 * libcairo -- optional, to write out the colour map
 * libsdl2 -- optional, to compile mflock with live view (`mflock_sdl`)
 * parallel
 * liblua5.4.5 (included in src/lua-5.3.5/ )

 On Ubuntu 22.04 these packages can be installed by

  ``` shell
  # Ubuntu
  sudo apt-get install libcairo-dev
  sudo apt-get install libreadline-dev
  sudo apt install pkg-config
  sudo apt-get install zlib1g-dev
  sudo apt-get install libsdl2-dev
  sudo apt-get install parallel
  ```

* lua, included but has to be built, on a Linux machine that corresponds to:
  ```
  cd src/lua-5.3.5/
  # All possible platforms: aix bsd c89 freebsd generic linux macosx mingw posix solaris
  # e.g., on linux use:
  make linux
  ```
  Please check the Lua documentation if you are on another platform.


## Custom install

If your system does not support deb files you can still install
chromflock. Typically by

``` shell
make
sudo make install
```

please check what the **makefile** does before doing this. If you
change the **DESTDIR** OR **PREFIX** variables please update the
script called **chromflock** to reflect these changes.
