#!/bin/bash

SCRIPTNAME="`basename $0`"

error()
{
    echo ""
    echo "   $SCRIPTNAME platform versionmin versionfull"
    echo "        e.g. % $SCRIPTNAME Linux64 2.7 2.7.2"
    echo ""
}

if test "x$1" = "x"; then
    error
    exit 1
fi

if test "x$2" = "x"; then
    error
    exit 1
fi

if test "x$3" = "x"; then
    error
    exit 1
fi

BUILD="$1"
PYVERSION="$2"
PYFULLVERSION="$3"

ARCH=`arch`
SRCPATH=`pwd`
DATE=`date '+%F-%T'`
SYS=`uname -s`

BUILDDIR="$SRCPATH/$BUILD/Python-$PYFULLVERSION/BUILD"
INSTALLDIR="$SRCPATH/$BUILD/Python-$PYFULLVERSION/INSTALL"
EGGDIR="$INSTALLDIR/eggs"

if [ ! -d $EGGDIR ]; then
    mkdir -p $EGGDIR
fi

PATH="$INSTALLDIR/bin:${PATH}"
LD_LIBRARY_PATH="$INSTALLDIR/lib"
PKG_CONFIG_PATH="$INSTALLDIR/lib/pkgconfig"
export LD_LIBRARY_PATH PATH PKG_CONFIG_PATH

if [ "$SYS" = "Darwin" ]; then
    TAR=gnutar
    CFLAGS="-isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5"
    RPATH="-Wl,-rpath,${PLACEHOLD}"

    DYLD_LIBRARY_PATH="$INSTALLDIR/lib"
    export DYLD_LIBRARY_PATH

    if [ "$BUILD" = "OSXI" ]
    then
	ARCHFLAGS="-arch i386"
	export ARCHFLAGS
    fi
fi

PYTHON="$INSTALLDIR/bin/python$PYVERSION"

TAR=tar
CC=gcc
CXX=g++
FC=gfortran


IPYTHONVER="0.11"
IPYTHON="ipython-$IPYTHONVER"
IPYTHONFILE="$SRCPATH/$IPYTHON.tar.gz"

NUMPYVER="1.6.1"
NUMPY="numpy-$NUMPYVER"
NUMPYFILE="$SRCPATH/$NUMPY.tar.gz"

MATPLOTLIBVER="1.1.0"
MATPLOTLIB="matplotlib-$MATPLOTLIBVER"
MATPLOTLIBFILE="$SRCPATH/$MATPLOTLIB.tar.gz"

PYFITSVER="3.0"
PYFITS="pyfits-$PYFITSVER"
PYFITSFILE="$SRCPATH/$PYFITS.tar.gz"

SETUPTOOLSVER="0.6c11"
SETUPTOOLS="setuptools-$SETUPTOOLSVER"
SETUPTOOLSFILE="$SRCPATH/$SETUPTOOLS.tar.gz"

DISTRIBUTEVER="0.6.21"
DISTRIBUTE="distribute-$DISTRIBUTEVER"
DISTRIBUTEFILE="$SRCPATH/$DISTRIBUTE.tar.gz"

ENSTALLERVER="4.4.1-1"
ENSTALLER="enstaller-$ENSTALLERVER"
ENSTALLERFILE="$SRCPATH/$ENSTALLER.tar.gz"

VOVER="0.7.2"
VO="vo-$VOVER"
VOFILE="$SRCPATH/$VO.tar.gz"

PYWCSVER="1.10-4.7"
PYWCS="pywcs-$PYWCSVER"
PYWCSFILE="$SRCPATH/$PYWCS.tar.gz"

ASCIITABLEVER="0.7.1"
ASCIITABLE="asciitable-$ASCIITABLEVER"
ASCIITABLEFILE="$SRCPATH/$ASCIITABLE.tar.gz"

PYPARSINGVER="1.5.6"
PYPARSING="pyparsing-$PYPARSINGVER"
PYPARSINGFILE="$SRCPATH/$PYPARSING.tar.gz"

PYREGIONVER="1.0.1"
PYREGION="pyregion-$PYREGIONVER"
PYREGIONFILE="$SRCPATH/$PYREGION.tar.gz"

ATPYVER="0.9.5.1"
ATPY="ATpy-$ATPYVER"
ATPYFILE="$SRCPATH/$ATPY.tar.gz"

APLPYVER="0.9.6"
APLPY="APLpy-$APLPYVER"
APLPYFILE="$SRCPATH/$APLPY.tar.gz"



make_eggsetup() {
    if [ ! -f setupegg.py ]; then
	echo \
	    "from setuptools import setup; " \
	    "execfile('setup.py')" \
	    > setupegg.py
    fi
}



build_setuptools()
{
    BUILDLOG="$BUILDDIR/build-$SETUPTOOLS-$DATE.log"

    cd $BUILDDIR
    $TAR -zxf $SETUPTOOLSFILE
    cd $SETUPTOOLS

    echo "Building $SETUPTOOLS ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $SETUPTOOLS"
        exit 1
    fi

    echo "Installing $SETUPTOOLS ..."
    $PYTHON setup.py install >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: installing $SETUPTOOLS"
        exit 1
    fi

    echo "Building $SETUPTOOLS egg ..."
    $PYTHON setup.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $SETUPTOOLS egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR
}


build_distribute()
{
    BUILDLOG="$BUILDDIR/build-$DISTRIBUTE-$DATE.log"

    cd $BUILDDIR
    $TAR -zxf $DISTRIBUTEFILE
    cd $DISTRIBUTE

    echo "Building $DISTRIBUTE ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $DISTRIBUTE"
        exit 1
    fi

    echo "Installing $DISTRIBUTE ..."
    $PYTHON setup.py install >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: installing $DISTRIBUTE"
        exit 1
    fi

    echo "Building $DISTRIBUTE egg ..."
    $PYTHON setup.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $DISTRIBUTE egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR
}


build_enstaller()
{
    BUILDLOG="$BUILDDIR/build-$ENSTALLER-$DATE.log"

    cd $BUILDDIR
    $TAR -zxf $ENSTALLERFILE
    cd $ENSTALLER

    echo "Building $ENSTALLER ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $ENSTALLER"
        exit 1
    fi

    echo "Building $ENSTALLER egg ..."
    $PYTHON setup.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $ENSTALLER egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR

    echo "Installing $ENSTALLER ..."
    $PYTHON setup.py install >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: installing $ENSTALLER"
        exit 1
    fi
}


build_ipython()
{
    BUILDLOG="$BUILDDIR/build-$IPYTHON-$DATE.log"

    cd $BUILDDIR
    $TAR -zxf $IPYTHONFILE
    cd $IPYTHON

    echo "Building $IPYTHON ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $IPYTHON"
        exit 1
    fi

    echo "Building $IPYTHON egg ..."
    $PYTHON setupegg.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $IPYTHON egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR

    # $PYTHON setup.py install >>$BUILDLOG 2>&1
    # if(( $? )); then
    #     echo "ERROR: installing $IPYTHON"
    #     exit 1
    # fi
}

build_numpy()
{
    BUILDLOG="$BUILDDIR/build-$NUMPY-$DATE.log"
    
    cd $BUILDDIR
    $TAR -zxf $NUMPYFILE
    cd $NUMPY

    echo "Building $NUMPY ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $NUMPY"
        exit 1
    fi

    echo "Building $NUMPY egg ..."
    $PYTHON setupegg.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $NUMPY egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR

    echo "Installing $NUMPY ..."
    $PYTHON setup.py install >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: installing $NUMPY"
        exit 1
    fi
}


build_pyfits()
{
    BUILDLOG="$BUILDDIR/build-$PYFITS-$DATE.log"
    
    cd $BUILDDIR
    $TAR -zxf $PYFITSFILE
    cd $PYFITS

    echo "Building $PYFITS ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $PYFITS"
        exit 1
    fi

    echo "Building $PYFITS egg ..."
    $PYTHON setup.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $PYFITS egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR
}




build_matplotlib()
{
    BUILDLOG="$BUILDDIR/build-$MATPLOTLIB-$DATE.log"
    
    cd $BUILDDIR
    $TAR -zxf $MATPLOTLIBFILE
    cd $MATPLOTLIB

    echo "Building $MATPLOTLIB ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $MATPLOTLIB"
        exit 1
    fi

    echo "Building $MATPLOTLIB egg ..."
    $PYTHON setupegg.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $MATPLOTLIB egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR
}


build_vo()
{
    BUILDLOG="$BUILDDIR/build-$VO-$DATE.log"
    
    cd $BUILDDIR
    $TAR -zxf $VOFILE
    cd $VO

    echo "Building $VO ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $VO"
        exit 1
    fi

    make_eggsetup
    
    echo "Building $VO egg ..."
    $PYTHON setupegg.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $VO egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR
}


build_pywcs()
{
    BUILDLOG="$BUILDDIR/build-$PYWCS-$DATE.log"
    
    cd $BUILDDIR
    $TAR -zxf $PYWCSFILE
    cd $PYWCS

    echo "Building $PYWCS ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $PYWCS"
        exit 1
    fi

    make_eggsetup

    echo "Building $PYWCS egg ..."
    $PYTHON setupegg.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $PYWCS egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR
}

build_asciitable()
{
    BUILDLOG="$BUILDDIR/build-$ASCIITABLE-$DATE.log"
    
    cd $BUILDDIR
    $TAR -zxf $ASCIITABLEFILE
    cd $ASCIITABLE

    echo "Building $ASCIITABLE ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $ASCIITABLE"
        exit 1
    fi

    make_eggsetup

    echo "Building $ASCIITABLE egg ..."
    $PYTHON setupegg.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $ASCIITABLE egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR
}


build_pyparsing()
{
    BUILDLOG="$BUILDDIR/build-$PYPARSING-$DATE.log"
    
    cd $BUILDDIR
    $TAR -zxf $PYPARSINGFILE
    cd $PYPARSING

    echo "Building $PYPARSING ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $PYPARSING"
        exit 1
    fi

    make_eggsetup

    echo "Building $PYPARSING egg ..."
    $PYTHON setupegg.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $PYPARSING egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR
}


build_pyregion()
{
    BUILDLOG="$BUILDDIR/build-$PYREGION-$DATE.log"
    
    cd $BUILDDIR
    $TAR -zxf $PYREGIONFILE
    cd $PYREGION

    echo "Building $PYREGION ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $PYREGION"
        exit 1
    fi

    make_eggsetup

    echo "Building $PYREGION egg ..."
    $PYTHON setup.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $PYREGION egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR
}


build_atpy()
{
    BUILDLOG="$BUILDDIR/build-$ATPY-$DATE.log"
    
    cd $BUILDDIR
    $TAR -zxf $ATPYFILE
    cd $ATPY

    echo "Building $ATPY ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $ATPY"
        exit 1
    fi
    
    make_eggsetup

    echo "Building $ATPY egg ..."
    $PYTHON setupegg.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $ATPY egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR
}


build_aplpy()
{
    BUILDLOG="$BUILDDIR/build-$APLPY-$DATE.log"
    
    cd $BUILDDIR
    $TAR -zxf $APLPYFILE
    cd $APLPY

    echo "Building $APLPY ..."
    $PYTHON setup.py build >$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $APLPY"
        exit 1
    fi

    make_eggsetup

    echo "Building $APLPY egg ..."
    $PYTHON setupegg.py bdist_egg >>$BUILDLOG 2>&1
    if(( $? )); then
        echo "ERROR: building $APLPY egg"
        exit 1
    fi

    cp dist/*.egg $EGGDIR
}


# SCRIPT

build_setuptools
build_distribute
build_enstaller
build_numpy
build_ipython
build_pyfits
build_matplotlib

# build_vo
# build_pywcs
# build_pyparsing
# build_pyregion
# build_asciitable
# build_atpy
# build_aplpy

exit 0
