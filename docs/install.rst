Installation
============

System requirements
-------------------

Only builds on x86-64 running Linux (possibly also BSD and MacOS).

For compilation the following is needed:

-  gcc or clang with OpenMP
-  make and cmake
-  libreadline – to compile lua
-  zlib – to write ``.gz`` files
-  pkg-config – for the makefile
-  libcairo – optional, to write out the colour map
-  libsdl2 – optional, enables the "live view" in mflock.
-  parallel - for distributing the mflock jobs over multiple cores
-  liblua5.3.5 (included in src/lua-5.3.5/ )

On Ubuntu 22.04 these packages can be installed by

.. code:: shell

   # Ubuntu
   sudo apt-get install libcairo-dev
   sudo apt-get install libreadline-dev
   sudo apt install pkg-config
   sudo apt-get install zlib1g-dev
   sudo apt-get install libsdl2-dev
   sudo apt-get install parallel

-  lua, included but has to be built, on a Linux machine that
   corresponds to:

   ::

      cd src/lua-5.3.5/
      # All possible platforms: aix bsd c89 freebsd generic linux macosx mingw posix solaris
      # e.g., on linux use:
      make linux

   Please check the Lua documentation if you are on another platform.

On Fedora 42

.. code:: shell

   sudo dnf install sdl2-compat-devel


Typical build and install
-------------------------

.. code:: shell

   cd src/lua-5.3.5/
   make linux # or pick another suitable target
   cd ../../
   mkdir build
   cd build
   cmake ..
   make
   sudo make install
