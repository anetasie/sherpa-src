#!/bin/bash

SCRIPTNAME="`basename $0`"

error() 
{
    echo ""
    echo "   $SCRIPTNAME system python"
    echo "        e.g. % $SCRIPTNAME Linux64 Python-2.7.2"
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
PYTHON="$2"
SYS=`uname -s`
SRC=`pwd`
DATE=`date '+%F-%T'`
TAR=tar

BUILDDIR="$SRC/$BUILD/$PYTHON/BUILD"
INSTALLDIR="$SRC/$BUILD/$PYTHON/INSTALL"

PLACEHOLD="/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD"

CFLAGS="-fPIC"
LDFLAGS="-fPIC"
ARCHFLAGS=
RPATH="-Wl,-R${PLACEHOLD}"

if [ $SYS = "Darwin" ]; then
    TAR=gnutar

    if [ "$BUILD" = "OSXI" ]; then
	ARCHFLAGS="-arch i386"
    fi

    CFLAGS="$CFLAGS -isysroot /Developer/SDKs/MacOSX10.6.sdk -mmacosx-version-min=10.6 $ARCHFLAGS"
    RPATH="-Wl,-rpath,${PLACEHOLD}"
    LDFLAGS="$CFLAGS"
    export ARCHFLAGS CFLAGS LDFLAGS
fi

# OTS 

CHRPATHFILE="$SRC/chrpath-0.13.tar.gz"
CHRPATHDIR="chrpath-0.13"

SQLITEDIR="sqlite-autoconf-3070701"
SQLITEFILE="$SRC/$SQLITEDIR.tar.gz"

TKFILE="$SRC/tk8.5.10-src.tar.gz"
TCLFILE="$SRC/tcl8.5.10-src.tar.gz"

TCLDIR="tcl8.5.10/unix"
TKDIR="tk8.5.10/unix"

READLINEFILE="$SRC/readline-6.2.tar.gz"
READLINEDIR="readline-6.2"

SSLFILE="$SRC/openssl-0.9.8e.tar.gz"
SSLDIR="openssl-0.9.8e"

ZLIBFILE="$SRC/zlib-1.2.5.tar.gz"
ZLIBDIR="zlib-1.2.5"

BZ2FILE="$SRC/bzip2-1.0.6.tar.gz"
BZ2DIR="bzip2-1.0.6"

NCURSESFILE="$SRC/ncurses-5.9.tar.gz"
NCURSESDIR="ncurses-5.9"

FREETYPEVER="2.4.6"
FREETYPE="freetype-$FREETYPEVER"
FREETYPEFILE="$SRC/$FREETYPE.tar.bz2"

PNGVER="1.2.46"
PNG="libpng-$PNGVER"
PNGFILE="$SRC/$PNG.tar.bz2"

DBVER="4.7.25"
DB="db-$DBVER"
DBFILE="$SRC/$DB.tar.gz"

PKGVER="0.23"
PKG="pkg-config-$PKGVER"
PKGFILE="$SRC/$PKG.tar.gz"


setup() 
{
    if [ ! -d $BUILDDIR ]; then
	mkdir -p $BUILDDIR
    fi

    if [ ! -d $INSTALLDIR ]; then
	mkdir -p $INSTALLDIR
    fi
}


make_chrpath()
{
    cd $BUILDDIR

    if [ ! -f $CHRPATHFILE ]; then
	echo "Error: $CHRPATHFILE not found!"
	exit 1
    fi

    $TAR -zxf $CHRPATHFILE
    if (( $? )); then
        echo "ERROR: Untarring $CHRPATHFILE failed!"
        exit 1
    fi    

    cd $CHRPATHDIR

    echo "Configuring $CHRPATHDIR ..."
    ./configure CFLAGS="$CFLAGS $ARCHFLAGS" LDFLAGS="$CFLAGS $ARCHFLAGS" \
	--prefix=$INSTALLDIR >/dev/null 2>&1
    if (( $? )); then
        echo "ERROR: configuring $CHRPATHDIR"
        exit 1
    fi

    echo "Building $CHRPATHDIR ..."
    make  >/dev/null 2>&1
    if (( $? )); then
        echo "ERROR: making $CHRPATHDIR"
        exit 1
    fi    

    echo "Installing $CHRPATHDIR ..."
    make install  >/dev/null 2>&1
    if (( $? )); then
        echo "ERROR: installing $CHRPATHDIR"
        exit 1
    fi
}

make_tcl()
{
    cd $BUILDDIR

    if [ ! -f $TCLFILE ]; then
	echo "Error: $TCLFILE not found!"
	exit 1
    fi

    $TAR -zxf $TCLFILE
    if (( $? )); then
        echo "ERROR: Untarring $TCLFILE failed!"
        exit 1
    fi 

    cd $TCLDIR

    echo "Configuring $TCLDIR ..."
    ./configure CFLAGS="$CFLAGS $ARCHFLAGS" LDFLAGS="$CFLAGS $ARCHFLAGS -L${INSTALLDIR}/lib ${RPATH}" \
	--prefix=$INSTALLDIR --enable-shared >"$BUILDDIR/tcl-config-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: configuring $TCLDIR"
        exit 1
    fi

    echo "Building $TCLDIR ..."
    make >"$BUILDDIR/tcl-build-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: making $TCLDIR"
        exit 1
    fi

    echo "Installing $TCLDIR ..."
    make install >"$BUILDDIR/tcl-install-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: installing $TCLDIR"
        exit 1
    fi
}


make_tk()
{
    cd $BUILDDIR

    if [ ! -f $TKFILE ]; then
	echo "Error: $TKFILE not found!"
	exit 1
    fi

    $TAR -zxf $TKFILE
    if (( $? )); then
        echo "ERROR: Untarring $TKFILE failed!"
        exit 1
    fi 

    cd $TKDIR

    echo "Configuring $TKDIR ..."
    ./configure CFLAGS="$CFLAGS $ARCHFLAGS" LDFLAGS="$CFLAGS $ARCHFLAGS -L${INSTALLDIR}/lib ${RPATH}" \
	--prefix=$INSTALLDIR --enable-shared --with-tcl=$BUILDDIR/$TCLDIR >"$BUILDDIR/tk-config-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: configuring $TKDIR"
        exit 1
    fi

    echo "Building $TKDIR ..."    
    make  >"$BUILDDIR/tk-build-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: making $TKDIR"
        exit 1
    fi
    
    echo "Installing $TKDIR ..."
    make install >"$BUILDDIR/tk-install-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: installing $TKDIR"
        exit 1
    fi
}

make_readline() 
{
    cd $BUILDDIR

    if [ ! -f $READLINEFILE ]; then
	echo "Error: $READLINEFILE not found!"
	exit 1
    fi

    $TAR -zxf $READLINEFILE
    if (( $? )); then
        echo "ERROR: Untarring $READLINEFILE failed!"
        exit 1
    fi 

    cd $READLINEDIR

    echo "Configuring $READLINEDIR ..."
    ./configure CFLAGS="$CFLAGS $ARCHFLAGS" LDFLAGS="$CFLAGS $ARCHFLAGS -L${INSTALLDIR}/lib ${RPATH}" \
	--prefix=$INSTALLDIR --enable-shared >"$BUILDDIR/readline-config-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: configuring $READLINEDIR"
        exit 1
    fi

    echo "Building $READLINEDIR ..."
    make >"$BUILDDIR/readline-build-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: making $READLINEDIR"
        exit 1
    fi

    echo "Installing $READLINEDIR ..."
    make install >"$BUILDDIR/readline-install-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: installing $READLINEDIR"
        exit 1
    fi

    if [ $SYS = "Darwin" ]; then
     	cd $INSTALLDIR/lib
	
     	/bin/rm -f libhistory*dylib*
     	/bin/rm -f libreadline*dylib*
    fi
}

make_sqlite()
{
    cd $BUILDDIR

    if [ ! -f $SQLITEFILE ]; then
	echo "Error: $SQLITEFILE not found!"
	exit 1
    fi

    $TAR -zxf $SQLITEFILE
    if (( $? )); then
        echo "ERROR: Untarring $SQLITEFILE failed!"
        exit 1
    fi 

    cd $SQLITEDIR

    echo "Configuring $SQLITEDIR ..."
    ./configure --prefix=$INSTALLDIR --enable-shared >"$BUILDDIR/sqlite3-config-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: configuring $SQLITEDIR"
        exit 1
    fi

    echo "Building $SQLITEDIR ..."
    make  >"$BUILDDIR/sqlite3-build-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: making $SQLITEDIR"
        exit 1
    fi
    

    echo "Installing $SQLITEDIR ..."
    make install  >"$BUILDDIR/sqlite3-install-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: installing $SQLITEDIR"
        exit 1
    fi
}

make_openssl()
{
    cd $BUILDDIR

    if [ ! -f $SSLFILE ]; then
	echo "Error: $SSLFILE not found!"
	exit 1
    fi

    $TAR -zxf $SSLFILE
    if (( $? )); then
        echo "ERROR: Untarring $SSLFILE failed!"
        exit 1
    fi 

    cd $SSLDIR

    SSLOPTIONS="shared no-asm no-sse2"
    if [ "$BUILD" = "Linux" ]; then
	SSLOPTIONS="$SSLOPTIONS 386"
    fi

    echo "Configuring $SSLDIR ..."    
    ./config --prefix=$INSTALLDIR --openssldir=$INSTALLDIR/openssl $SSLOPTIONS >"$BUILDDIR/openssl-config-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: configuring $SSLDIR"
        exit 1
    fi

    echo "Building $SSLDIR ..."
    make >"$BUILDDIR/openssl-build-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: making $SSLDIR"
        exit 1
    fi
    
    echo "Testing $SSLDIR ..."
    make test >"$BUILDDIR/openssl-test-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: testing $SSLDIR"
        exit 1
    fi

    echo "Installing $SSLDIR ..."
    make install >"$BUILDDIR/openssl-install-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: installing $SSLDIR"
        exit 1
    fi
}

make_zlib()
{
    # Zlib is built statically, linux only

    cd $BUILDDIR

    if [ ! -f $ZLIBFILE ]; then
	echo "Error: $ZLIBFILE not found!"
	exit 1
    fi

    $TAR -zxf $ZLIBFILE
    if (( $? )); then
        echo "ERROR: Untarring $ZLIBFILE failed!"
        exit 1
    fi 

    cd $ZLIBDIR

    echo "Configuring $ZLIBDIR ..."
    ./configure --prefix=$INSTALLDIR >"$BUILDDIR/zlib-config-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: configuring $ZLIBDIR"
        exit 1
    fi

    echo "Building $ZLIBDIR ..."
    make >"$BUILDDIR/zlib-build-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: making $ZLIBDIR"
        exit 1
    fi
    

    echo "Installing $ZLIBDIR ..."
    make install >"$BUILDDIR/zlib-install-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: installing $ZLIBDIR"
        exit 1
    fi
}

make_bzip2()
{
    # Bzip2 is built statically, linux only

    cd $BUILDDIR

    if [ ! -f $BZ2FILE ]; then
	echo "Error: $BZ2FILE not found!"
	exit 1
    fi

    $TAR -zxf $BZ2FILE
    if (( $? )); then
        echo "ERROR: Untarring $BZ2FILE failed!"
        exit 1
    fi 

    cd $BZ2DIR

    echo "Building $BZ2DIR ..."
    make -f Makefile-libbz2_so >"$BUILDDIR/bzip2-build-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: making $BZ2DIR"
        exit 1
    fi

    echo "Installing $BZ2DIR ..."
    make PREFIX=$INSTALLDIR install >"$BUILDDIR/bzip2-install-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: installing $BZ2DIR"
        exit 1
    fi
}

make_ncurses()
{
    cd $BUILDDIR

    if [ ! -f $NCURSESFILE ]; then
	echo "Error: $NCURSESFILE not found!"
	exit 1
    fi

    $TAR -zxf $NCURSESFILE
    if (( $? )); then
        echo "ERROR: Untarring $NCURSESFILE failed!"
        exit 1
    fi 

    cd $NCURSESDIR

    echo "Configuring $NCURSESDIR ..."
    ./configure CFLAGS="$CFLAGS $ARCHFLAGS" LDFLAGS="$CFLAGS $ARCHFLAGS -L${INSTALLDIR}/lib ${RPATH}" \
	--prefix=$INSTALLDIR \
	--without-normal \
        --enable-shared \
	--without-cxx \
        --with-shared \
	--enable-termcap \
	--without-debug >"$BUILDDIR/ncurses-config-$DATE.log" 2>&1

    if (( $? )); then
        echo "ERROR: configuring $NCURSESDIR"
        exit 1
    fi

    echo "Building $NCURSESDIR ..."
    make  >"$BUILDDIR/ncurses-build-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: making $NCURSESDIR"
        exit 1
    fi    

    echo "Installing $NCURSESDIR ..."
    make install >"$BUILDDIR/ncurses-install-$DATE.log" 2>&1
    if (( $? )); then
        echo "ERROR: installing $NCURSESDIR"
        exit 1
    fi
}

make_freetype()
{
   
    cd $BUILDDIR

    if [ ! -f $FREETYPEFILE ]; then
	echo "Error: $FREETYPEFILE not found!"
	exit 1
    fi

    $TAR -jxf $FREETYPEFILE
    if (( $? )); then
        echo "ERROR: Untarring $FREETYPEFILE failed!"
        exit 1
    fi 

    cd $FREETYPE

    echo "Configuring $FREETYPE ..."
    ./configure CFLAGS="$CFLAGS $ARCHFLAGS" LDFLAGS="$CFLAGS $ARCHFLAGS -L${INSTALLDIR}/lib ${RPATH}" \
	--prefix=$INSTALLDIR >"$BUILDDIR/freetype-config-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: configuring $FREETYPE"
        exit 1
    fi

    echo "Building $FREETYPE ..."
    make >"$BUILDDIR/freetype-build-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: making $FREETYPE"
        exit 1
    fi

    echo "Installing $FREETYPE ..."
    make install >"$BUILDDIR/freetype-install-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: installing $FREETYPE"
        exit 1
    fi
}


make_libpng()
{   
    cd $BUILDDIR

    if [ ! -f $PNGFILE ]; then
	echo "Error: $PNGFILE not found!"
	exit 1
    fi

    $TAR -jxf $PNGFILE
    if (( $? )); then
        echo "ERROR: Untarring $PNGFILE failed!"
        exit 1
    fi 

    cd $PNG

    echo "Configuring $PNG ..."
    ./configure CFLAGS="$CFLAGS $ARCHFLAGS" LDFLAGS="$CFLAGS $ARCHFLAGS -L${INSTALLDIR}/lib ${RPATH}" \
	--prefix=$INSTALLDIR >"$BUILDDIR/libpng-config-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: configuring $PNG"
        exit 1
    fi

    echo "Building $PNG ..."
    make >"$BUILDDIR/libpng-build-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: making $PNG"
        exit 1
    fi

    echo "Installing $PNG ..."
    make install >"$BUILDDIR/libpng-install-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: installing $PNG"
        exit 1
    fi
}


make_db()
{   
    cd $BUILDDIR
    if [ ! -f $DBFILE ]; then
	echo "Error: $DBFILE not found!"
	exit 1
    fi

    $TAR -zxf $DBFILE
    if (( $? )); then
        echo "ERROR: Untarring $DBFILE failed!"
        exit 1
    fi 
    cd $DB/build_unix

    echo "Configuring $DB ..."
    ../dist/configure CFLAGS="$CFLAGS $ARCHFLAGS" LDFLAGS="$CFLAGS $ARCHFLAGS -L${INSTALLDIR}/lib ${RPATH}" --includedir=$INSTALLDIR/include/db4 \
	--prefix=$INSTALLDIR --disable-shared >"$BUILDDIR/db-config-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: configuring $DB"
        exit 1
    fi

    echo "Building $DB ..."
    make >"$BUILDDIR/db-build-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: making $DB"
        exit 1
    fi

    echo "Installing $DB ..."
    make install >"$BUILDDIR/db-install-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: installing $DB"
        exit 1
    fi


}


make_pkgconfig()
{

    cd $BUILDDIR
    if [ ! -f $PKGFILE ]; then
	echo "Error: $PKGFILE not found!"
	exit 1
    fi

    $TAR -zxf $PKGFILE
    if (( $? )); then
        echo "ERROR: Untarring $PKGFILE failed!"
        exit 1
    fi 
    cd $PKG

    echo "Configuring $PKG..."
    ./configure CFLAGS="$CFLAGS $ARCHFLAGS" LDFLAGS="$CFLAGS $ARCHFLAGS -L${INSTALLDIR}/lib ${RPATH}" --prefix=$INSTALLDIR >"$BUILDDIR/pkg-config-config-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: configuring $PKG"
        exit 1
    fi

    echo "Building $PKG ..."
    make >"$BUILDDIR/pkg-config-build-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: making $PKG"
        exit 1
    fi

    echo "Installing $PKG ..."
    make install >"$BUILDDIR/pkg-config-install-$DATE.log" 2>&1
    if(( $? )); then
        echo "ERROR: installing $PKG"
        exit 1
    fi
}


# Script

if [ -f $INSTALLDIR/bin/python ]; then  
    echo "Python found! Skipping OTS build..."
    exit 0
fi

setup

if [ ! -f $INSTALLDIR/lib/libdb.a ]; then
    make_db
fi

# rely on MacOS for framework Tk and Tcl
if [ \( $BUILD = "Linux" \) -o \( $BUILD = "Linux64" \) ]; then

    if [ ! -f $INSTALLDIR/lib/tclConfig.sh ]; then
	make_tcl
    fi
    
    if [ ! -f $INSTALLDIR/lib/tkConfig.sh ]; then
	make_tk
    fi
fi

if [ ! -f $INSTALLDIR/lib/libreadline.a ]; then
    make_readline
fi

if [ \( $BUILD = "Linux" \) -o \( $BUILD = "Linux64" \) ]; then

    if [ ! -f $INSTALLDIR/bin/chrpath ]; then
	make_chrpath
    fi

    if [ ! -f $INSTALLDIR/lib/libz.a ]; then
	make_zlib
    fi

    if [ ! -f $INSTALLDIR/lib/libbz2.a ]; then
	make_bzip2
    fi

    if [ ! -f $INSTALLDIR/lib/libsqlite3.a ]; then
	make_sqlite
    fi

    if [ ! -f $INSTALLDIR/lib/libssl.a ]; then
	make_openssl
    fi

    if [ ! -f $INSTALLDIR/lib/libncurses.so ]; then
	make_ncurses
    fi
fi

if [ ! -f $INSTALLDIR/lib/libfreetype.a ]; then
    make_freetype
fi

# EPD relies on system libpng12 on Linux
# EPD builds libpng12 AND libjpeg on Mac

# CentOS does not have libpng by default

#if [ \( $BUILD = "OSXI" \) -o \( $BUILD = "OSX64" \) ]; then
#fi

if [ ! -f $INSTALLDIR/lib/libpng12.a ]; then
    make_libpng
fi

if [ ! -f $INSTALLDIR/bin/pkg-config ]; then
    make_pkgconfig
fi

exit 0
