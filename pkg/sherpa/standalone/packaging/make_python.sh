#!/bin/bash

SCRIPTNAME=`basename $0`
SCRIPTDIR=`dirname $0`
SCRIPTDIR=`(cd $SCRIPTDIR; pwd)`

error()
{
    echo ""
    echo "   $SCRIPTNAME platform version"
    echo "        e.g. % $SCRIPTNAME Linux64 2.7.2"
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

BUILD="$1"
VERSION="$2"
MAJOR=`echo $VERSION | awk '{print substr($1,1,3)}'`
ARCH=`arch`
SRCPATH=`pwd`
DATE=`date '+%F-%T'`
SYS=`uname -s`

PYTHON="Python-$VERSION"
PYTHONFILE="$SRCPATH/$PYTHON.tgz"

TAR=tar
CC=gcc

BUILDDIR="$SRCPATH/$BUILD/$PYTHON/BUILD"
INSTALLDIR="$SRCPATH/$BUILD/$PYTHON/INSTALL"

CONFIGLOG="$BUILDDIR/config-$PYTHON-$BUILD-$DATE.log"
BUILDLOG="$BUILDDIR/build-$PYTHON-$BUILD-$DATE.log"
INSTALLLOG="$BUILDDIR/install-$PYTHON-$BUILD-$DATE.log"

PLACEHOLD="/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD"


CFLAGS=
ARCHFLAGS=

RPATH="-Wl,-R${PLACEHOLD}"

LD_LIBRARY_PATH="$INSTALLDIR/lib"
export LD_LIBRARY_PATH

if [ "$SYS" = "Darwin" ]; then
    TAR=gnutar
    CFLAGS="-isysroot /Developer/SDKs/MacOSX10.6.sdk -mmacosx-version-min=10.6"
    RPATH="-Wl,-rpath,${PLACEHOLD}"

    DYLD_LIBRARY_PATH="$INSTALLDIR/lib"
    export DYLD_LIBRARY_PATH

    if [ "$BUILD" = "OSXI" ]
    then
	ARCHFLAGS="-arch i386"
	export ARCHFLAGS
    fi
fi

build_python()
{

    cd $BUILDDIR
    if [ ! -f $PYTHONFILE ]; then
	echo "ERROR: $PYTHONFILE not found!"
	exit 1
    fi

    echo "Untarring $PYTHON ..."
    $TAR -zxf $PYTHONFILE
    if (( $? )); then
        echo "ERROR: unzipping $PYTHONFILE"
        exit 1
    fi

    echo "Configuring $PYTHON ..."
    cd $PYTHON
    ./configure \
	OPT="-I${INSTALLDIR}/include" \
	LDFLAGS="-L${INSTALLDIR}/lib ${RPATH}" \
	--prefix=$INSTALLDIR \
	--enable-shared > $CONFIGLOG 2>&1

    if (( $? )); then
        echo "ERROR: configuring $PYTHON"
        exit 1
    fi

    echo "Making $PYTHON ..."
    make > $BUILDLOG 2>&1
    if (( $? )); then
        echo "ERROR: make $PYTHON"
        exit 1
    fi

    cd $SRCPATH
}


build_python_darwin()
{

    cd $BUILDDIR
    if [ ! -f $PYTHONFILE ]; then
	echo "ERROR: $PYTHONFILE not found!"
	exit 1
    fi

    echo "Untarring $PYTHON ..."
    $TAR -zxf $PYTHONFILE
    if (( $? )); then
        echo "ERROR: unzipping $PYTHONFILE"
        exit 1
    fi
    cd $PYTHON
 
    # edit setup.py to allow linking to the tk and tcl dylibs (Python expects tk and tcl frameworks)
    if [ "$BUILD" = "OSX64" ]; then
	echo "Editing setup.py to link tk and tcl as dylibs ..."
	mv setup.py setup.py.old
	sed 's/if (platform == '\''darwin'\'' and/if (False and/g' <setup.py.old >setup.py
    fi

    #echo "Editing Makefile.pre.in to link readline as dylib ..."
    # edit Makefile.pre.in to include a target for readline and termcap?
    #mv Makefile.pre.in Makefile.pre.in.old
    #sed 's/-current_version \$(VERSION);/-current_version \$(VERSION) \$(OTHER_LIBTOOL_OPT);/g' <Makefile.pre.in.old >Makefile.pre.in


    #--enable-framework=${INSTALLDIR}

    echo "Configuring $PYTHON ..."
    ./configure \
	OPT="-g -O3 -I${INSTALLDIR}/include" \
	CFLAGS="$CFLAGS $ARCHFLAGS" \
	LDFLAGS="$CFLAGS $ARCHFLAGS -L${INSTALLDIR}/lib ${RPATH}" \
	--enable-shared \
	--prefix=$INSTALLDIR > $CONFIGLOG 2>&1

    if (( $? )); then
        echo "ERROR: configuring $PYTHON"
        exit 1
    fi

    # edit Modules/Setup.local to include readline
    #echo readline readline.c -I${INSTALLDIR}/include -L${INSTALLDIR}/lib \
    #	-lreadline -ltermcap >> Modules/Setup.local

    echo "Making $PYTHON ..."
    #make OTHER_LIBTOOL_OPT="-L${INSTALLDIR}/lib -lreadline -lncurses" \
#	>$BUILDLOG 2>&1
    make >$BUILDLOG 2>&1
    if (( $? )); then
        echo "ERROR: make $PYTHON"
        exit 1
    fi

    cd $SRCPATH
}


install_python()
{
    cd $BUILDDIR

    cd $PYTHON
    echo "Installing $PYTHON ..."
    make install > $INSTALLLOG 2>&1
    if (( $? )); then
        echo "ERROR: make install $PYTHON"
        exit 1
    fi

    cd $SRCPATH
}


install_python_darwin()
{
    cd $BUILDDIR

    cd $PYTHON
    echo "Installing $PYTHON ..."
    #make install DESTDIR=${INSTALLDIR} >$INSTALLLOG 2>&1
    make install PYTHONAPPSDIR=${INSTALLDIR} >$INSTALLLOG 2>&1
    if (( $? )); then
        echo "ERROR: make install $PYTHON"
        exit 1
    fi

    cd $SRCPATH

    # edit binaries using install_name_tool

    PYLIB="libpython$MAJOR.dylib"
    PYLIBDIR="$INSTALLDIR/lib/$PYLIB"
    install_name_tool -change $PYLIBDIR $PYLIB $INSTALLDIR/bin/python
    install_name_tool -change $PYLIBDIR $PYLIB $INSTALLDIR/bin/python$MAJOR
    #install_name_tool -change $PYLIBDIR $PYLIB $ROOT/bin/python$PYVERSION-config
    chmod +w $PYLIBDIR
    install_name_tool -id $PYLIB $PYLIBDIR

}


# SCRIPT

if [ ! -d $BUILDDIR ]; then
    mkdir -p $BUILDDIR
fi

if [ ! -d $INSTALLDIR ]; then
    mkdir -p $INSTALLDIR
fi

if [ ! -f $BUILDDIR/$PYTHON/python ]; then
    if [ "$SYS" = "Darwin" ]; then
	build_python_darwin
    else
	build_python
    fi
fi

if [ ! -f $INSTALLDIR/bin/python ]; then
    if [ "$SYS" = "Darwin" ]; then
	install_python_darwin
    else
	install_python
    fi
fi

unset LD_LIBRARY_PATH DYLD_LIBRARY_PATH

exit 0
