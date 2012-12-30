Building pi-qmc
===============

Quick start
-----------

The easiest way to build is to use:

.. code-block:: bash

   ./configure make

For a parallel build

.. code-block:: bash

  ./configure --enable-mpi MPICXX=mpic++ MPICC=mpicc MPIF77=mpif77

where you should use the names of your MPI enabled compilers.

You can also build for different numbers of physical dimensions (default is NDIM=3)

.. code-block:: bash

  ./configure --with-ndim=2

Required libraries
------------------

We use the following libraries in the pi code:

*   `libxml2`_

*   `blitz++`_

*   `hdf5`_

*   `fftw3`_

*   `BLAS`_ / `LAPACK`_

*   `gsl`_

Advanced build using multiple directories
-----------------------------------------

In research, we often want different versions of the executables, for example, versions with and
without MPI, or versions compiled for two-dimensional systems. To accomplish this, we make a
``pibuilds`` directory beside our svn checkout directory (``pi`` or ``pi-qmc``). We then make empty
subdirectories for each build, for example ``ndim2mpi`` for a two dimensional MPI version. A typical
directory structure is:

::

    codes/
      pi-qmc/
        configure
        src/
        lib/
      pibuilds/
        ndim1/
        ndim2/
        ndim3/
        ndim1mpi/
        ndim2mpi/
        ndim3mpi/
        debug/

..


To build, go into the empty build directory,

.. code-block:: bash

   cd ~/codes/pibuilds/ndim2mpi

Then run the configure script with the desired options

.. code-block:: bash

   ../../configure --with-ndim=2 --enable-mpi

You will probably want more configure options; see the platform specific instructions below for some
examples.

Then, make the code in that directory,

.. code-block:: bash

   make -j2

For conveniance, you can make a soft link to the executable

.. code-block:: bash

    ln -sf ~/codes/pibuilds/ndim3mpi ~/bin/pi2Dmpi

Platform specific instructions
------------------------------

Mac OS X
````````

All the dependencies are available through [http://www.macports.org/ macports]. It is also handy to
install the latest gcc compilers (with gfortran), openmpi, and python utilities for data analysis and
plotting.

::

    $ port installed
      libxml2 @2.7.3_0 (active)
      blitz @0.9_0 (active)
      hdf5-18 @1.8.3_0 (active)
      gsl @1.12_0 (active)
      fftw-3 @3.2.2_0 (active)
    
      gcc44 @4.4.0_0 (active)
    
      python26 @2.6.2_3 (active)
      py26-numpy @1.3.0_0 (active)
      py26-ipython @0.9.1_0+scientific (active)
      py26-scipy @0.7.0_0+gcc44 (active)
      py26-tables @2.1_0 (active)

..

A bash script to download the pi-qmc source from git hub, compile on Mac OS X,
and run all tests is included in ``doc/deploy/macosx/start.sh``:

.. literalinclude:: ../../doc/deploy/macosx/start.sh
   :language: bash


The following configure works well on an intel mac:

::

    ../../pi/configure CXX=g++-mp-4.7 CC=gcc-mp-4.7 \
      CXXFLAGS="-O3 -g -Wall -ffast-math -ftree-vectorize \
      -march=native -fomit-frame-pointer -pipe" \ 
      F77=gfortran-mp-4.7

..



or, for an MPI enabled build,

::

    ../../pi/configure --enable-mpi CXX=g++-mp-4.4 CC=gcc-mp-4.4 F77=gfortran-mp-4.4 \
      MPICC=openmpicc MPICXX=openmpicxx MPIF77=openmpif77 \
      CXXFLAGS="-O3 -g -Wall -ffast-math -ftree-vectorize \
      -march=native -fomit-frame-pointer -pipe"

..



On a G5 mac, try:

::

    ../../pi/configure --with-ndim=3  F77=gfortran-mp-4.4 CC=gcc-mp-4.4 CXX=g++-mp-4.4\
      CXXFLAGS="-g -O3 -ffast-math -ftree-vectorize -maltivec -mpowerpc-gpopt \
      -mpowerpc64 falign-functions=32 -falign-labels=32 -falign-loops=32 -falign-jumps=32 -funroll-loops"

..



or, for an MPI enabled build,

::

    ../../pi/configure --with-ndim=3 --enable-mpi \
      CXXFLAGS="-g -O3 -ffast-math -ftree-vectorize -maltivec -mpowerpc-gpopt \
      -mpowerpc64 falign-functions=32 -falign-labels=32 -falign-loops=32 -falign-jumps=32 -funroll-loops" \
      F77=gfortran-mp-4.4 CC=gcc-mp-4.4 CXX=g++-mp-4.4  MPICC=openmpicc MPICXX=openmpicxx MPIF77=openmpif77

..



Linux (CentOS 5.3)
``````````````````

You can download dependencies using ``yum``. First, you may need to add access to the fedora
[http://fedoraproject.org/wiki/EPEL Extra Packages for Enterprise Linux (EPEL)].

sudo rpm -Uvh http://download.fedora.redhat.com/pub/epel/5/i386/epel-release-5-3.noarch.rpm

Then install the required packages for *_pi_*. (You probably want to compile atlas yourself to get
automatic performance tuning for your hardware, but the yum install will work if you are impatient.)
Note: replace ``x86_64`` with ``i386`` if you are on a 32 bit machine.

::

    sudo yum install libxml2-devel-versionXXX.x86_64 (here I don't know the correct version)
    sudo yum install blitz-devel.x86_64
    sudo yum install fftw3-devel.x86_64
    sudo yum install hdf5-devel.x86_64
    sudo yum install atlas-sse3-devel.x86_64
    sudo yum install lapack-devel.x86_64
    sudo yum install gsl-devel.x86_64

..



It is useful to install the gcc 4.3 compilers.

::

    sudo yum install gcc43.x86_64
    sudo yum install gcc43-c++.x86_64
    sudo yum install gcc43-gfortran.x86_64

..



Also, you will want an MPI implementation if you want to run in parallel,

sudo yum install openmpi-devel.x86_64

The openmpi package will require that you run ``mpi-selector`` and open a new terminal to get the
executables. Use the ``mpi-selector --list`` option to see what is available, then set a system-wide
default.

sudo mpi-selector --system --set openmpi-1.2.7-gcc-x86_64

When you configure pi, you will probably need to specify the location of your BLAS and LAPACK routines,

::

    ../../pi/configure CXX=g++43 CC=gcc43 F77=gfortran43 CXXFLAGS=\
    "-g -O3 -ffast-math -ftree-vectorize -march=native -fomit-frame-pointer -pipe"\
     --with-blas="-L/usr/lib64/atlas -llapack -lf77blas"

..



For mpi, just add --enable-mpi.

For the python analysis utilities, you'll want to install ipython and matplotlib.

::

    sudo yum install python-matplotlib
    sudo yum install ipython
    sudo yum install scipy

..


The python package `pytables`_ for reading HDF5 files is also required for the analysis scripts, but it
is not available through yum, so you'll have to download it and install it yourself.

HPC Centers
-----------

ASU Fulton: saguaro
```````````````````

For a serial build in two dimensions,

../../pi/configure --with-ndim=2 --enable-sprng CXX=icpc CC=icc CXXFLAGS="-O3 -xP -ipo" \
--with-blas="-L$MKL_LIB -lmkl_lapack -lmkl_intel_lp64 -lmkl_sequential -lmkl_core" \ F77=ifort AR="xild
-lib"

or for a parallel version,

../../pi/configure --with-ndim=2 --enable-sprng --enable-mpi MPICC=mpicc MPICXX=mpicxx \ CXX=icpc
CC=icc F77=ifort CXXFLAGS="-O3 -xP -ipo" AR="xild -lib" \ --with-blas="-L$MKL_LIB -lmkl_lapack
-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"

Omit the --enable-sprng option if you do not have the SPRNG library.

LONI-LSU: queenbee
``````````````````

You need to add some lines to your ``.soft`` file to include some required libraries,

::

    #My additions (CPATH mimics -I include directories).
    CPATH += /usr/local/packages/hdf5-1.8.1-intel10.1/include
    +gsl-1.9-intel10.1
    +sprng4-mvapich-1.1-intel-10.1
    +fftw-3.1.2-intel10.1
    CPATH += :/usr/local/packages/fftw-3.1.2-intel10.1/include
    +intel-mkl
    CPPFLAGS += -DMPICH_IGNORE_CXX_SEEK

..



For an MPI build, use,

::

    ../../pi/configure --with-ndim=3 --enable-mpi MPICC=mpicc MPICXX=mpicxx \
      CXX=icpc CC=icc F77=ifort AR="xild -lib" CXXFLAGS="-O3 -xP -ipo" \
      --with-blas="-lmkl_lapack -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"

..



NCSA: abe
`````````

You need to add some lines to your ``.soft`` file to include some required libraries,

::

    #My additions (CPATH mimics -I include directories).
    +libxml2-2.6.29
    +libxml2
    +intel-mkl
    +gsl-intel
    +hdf5-1.8.2
    CPATH += :/usr/apps/hdf/hdf5/v182/include
    LD_LIBRARY_PATH += /usr/apps/hdf/szip/lib
    +fftw-3.1-intel
    LD_LIBRARY_PATH += /usr/apps/math/fftw/fftw-3.1.2/intel10/lib
    CPATH += :/usr/apps/math/fftw/fftw-3.1.2/intel10/include
    +intel-mkl
    CPPFLAGS = "${CPPFLAGS} -DMPICH_IGNORE_CXX_SEEK"
    Also have blitz installed locally with --prefix=(your dir choice)

..



For an MPI build, use,

../../pi/configure --with-ndim=3 --enable-mpi MPICC=mpicc MPICXX=mpicxx CXX=icpc CC=icc \ CXXFLAGS="-O3
-xP -ipo" LDFLAGS="-lsz" \ --with-blas="-lmkl_lapack -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
F77=ifort AR="xild -lib"

TACC: Ranger
````````````

Cornell CNF: nanolab
````````````````````

The svn client wasn't working for me, so I built one in my ~/packages/bin directory. You need to
specify the most recent C++ and Fortran compilers by including the following in your .bash_profile,

::

    # Version 10 compilers
    source /opt/intel/cc/10.1.017/bin/iccvars.sh
    source /opt/intel/fc/10.1.017/bin/ifortvars.shsource /opt/intel/idb/10.1.017/bin/idbvars.sh
    source /opt/intel/mkl/10.0.4.023/tools/environment/mklvars32.sh

..



Also, make sure that ``/usr/lam-7.4.1_intelv10/bin`` is in your path to get the correct MPI compilers.

You need to build blitz (again, in my ~/packages directory). For a serial pi build,

::

    ../../pi/configure --with-ndim=3 CXX=icpc CC=icc CXXFLAGS="-O3 -ipo" \
    --with-blas="-Wl,-rpath,$MKLROOT/lib/32 -L/opt$MKLROOT/lib/32 -lmkl_intel \
    -lmkl_sequential -lmkl_core -lpthread -lm" F77=ifort AR="xild -lib"

..



::

    ../../pi/configure --with-ndim=3 CXX=icpc CC=icc CXXFLAGS="-O3 -ipo" \
    --with-blas="-Wl,-rpath,$MKLROOT/lib/32 -L$MKLROOT/lib/32 -lmkl_intel \
    -lmkl_sequential -lmkl_core -lpthread -lm" F77=ifort AR="xild -lib" \
    --enable-mpi MPICXX=mpic++ MPICC=mpicc MPIF77=mpif77

..



.. _`pytables`:
    http://www.pytables.org/

.. _`hdf5`:
    http://www.hdfgroup.org/

.. _`libxml2`:
    http://xmlsoft.org/

.. _`fftw3`:
    http://www.fftw.org/

.. _`blitz++`:
    http://sourceforge.net/projects/blitz/

.. _`gsl`:
    http://www.gnu.org/software/gsl/

.. _`BLAS`:
    http://www.netlib.org/blas/

.. _`LAPACK`:
    http://www.netlib.org/lapack/


.. _`contents`:
    index.xhtml
