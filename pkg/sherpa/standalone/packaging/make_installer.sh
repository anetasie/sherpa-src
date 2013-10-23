#!/bin/bash

SCRIPTNAME="`basename $0`"

error()
{
    echo ""
    echo "   $SCRIPTNAME root name version fullversion path_to_install.sh"
    echo "        e.g. % $SCRIPTNAME INSTALL Python-2.7.2 2.7 2.7.2 ../../install.sh"
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

if test "x$4" = "x"; then
    error
    exit 1
fi

if test "x$5" = "x"; then
    error
    exit 1
fi

SRC=`pwd`
DIR="$1"
NAME="$2"
NAMEFILE="../../../$NAME.sh"
PYVERSION="$3"
PYFULLVERSION="$4"
INSTALLSCRIPT="$5"



if [ ! -d $SRC/$DIR ]; then
    echo "Error: no such directory $DIR..."
    exit 1
fi

if [ ! -f $INSTALLSCRIPT ]; then
    echo "Error: $INSTALLSCRIPT does not exist..."
    exit 1
fi

ARCH=`arch`
SRCPATH=`pwd`
DATE=`date '+%F-%T'`
SYS=`uname -s`

TAR=tar
if [ "$SYS" = "Darwin" ]; then
    TAR=gnutar
fi

TARFILE="python.tar"
PAYLOADFILE="payload.tar"
PAYLOADFILEGZ="payload.tar.gz"

build_installer()
{

    cd $DIR
    
    mv ../$TARFILE .

    ARGS="$TARFILE"
    if [ -d scripts ]; then
	ARGS="$ARGS scripts"
    fi

    if [ -d eggs ]; then
	ARGS="$ARGS eggs"
    fi

    $TAR -cf $PAYLOADFILE $ARGS
    if [ $? -ne 0 ]; then
	echo "Error: Tarring $NAME payload failed!"
	exit 1
    fi
 
    gzip $PAYLOADFILE
    if [ $? -ne 0 ]; then
	echo "Error: Gzipping $NAME payload failed!"
	exit 1
    fi
    
    sed "s/PYFULLVERSION/${PYFULLVERSION}/g" <$INSTALLSCRIPT >tmp.1
    sed "s/PYVERSION/${PYVERSION}/g" <tmp.1 >tmp.2
    cat tmp.2 $PAYLOADFILEGZ >$NAMEFILE

    chmod +x $NAMEFILE

    /bin/rm -f $PAYLOADFILEGZ tmp.1 tmp.2
    echo "Created $NAME.sh"
}

# SCRIPT

if [ ! -f $NAMEFILE ]; then
    build_installer
fi

exit 0
