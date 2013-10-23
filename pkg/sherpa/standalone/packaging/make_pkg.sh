#!/bin/bash

SCRIPTNAME="`basename $0`"

if test "x$1" = "x"; then
    echo "$SCRIPTNAME root [ name version installdir ]"
    exit 1
fi

DIR="$1"
NAME="$DIR"
VERSION="$DIR"
INSTALL="$DIR"

if test "x$2" != "x"; then
    NAME="$2"
fi

if test "x$3" != "x"; then
    VERSION="$3"
fi

if test "x$4" != "x"; then
    INSTALL="$4"
fi

PKG="$NAME.mpkg"
SYS=`uname -s`
DATE=`date '+%F-%T'`

PKGLOG="pkg-$NAME-$DATE.log"

if [ "$SYS" != "Darwin" ]; then
    echo "Error: System must be Mac OS"
    exit 2
fi

if [ ! -d scripts ]; then
    mkdir -p scripts
fi

if [ ! -d resources ]; then
    mkdir -p resources
fi


build_pkg()
{
    /Developer/usr/bin/packagemaker \
	--title "$NAME installer" \
	--version "$VERSION" \
	--install-to "$INSTALL" \
	--resources resources/ \
	--scripts scripts/ \
	--root-volume-only \
	--domain system \
	--verbose \
	--target 10.5 \
	--id "edu.harvard.cxc.$NAME" \
	--root "$DIR" \
	--out "$PKG" \
	>$PKGLOG 2>&1

    if [ $? -ne 0 ] ; then
	echo "Error: Building package $NAME failed..."
	exit 3
    fi
}

echo "Creating $NAME MacOSX package..."
build_pkg 

exit 0
