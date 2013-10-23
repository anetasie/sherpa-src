#!/bin/bash

SCRIPTNAME="`basename $0`"

error() 
{
    echo "$SCRIPTNAME pkgfile dmgname"
}


if test "x$1" = "x"; then
    error
    exit 1
fi

if test "x$2" = "x"; then
    error
    exit 1
fi

PKG="$1"
NAME="$2"
SYS=`uname -s`
DATE=`date '+%F-%T'`

DMGLOG="dmg-$NAME-$DATE.log"

if [ "$SYS" != "Darwin" ]; then
    echo "Error: System must be Mac OS"
    exit 2
fi

build_dmg()
{

    DMG="$NAME-$DATE.dmg"
    OUT="$NAME.dmg"
    VOL="$NAME"
    
    echo "Creating $DMG of $PKG ..."
    hdiutil create \
	-format UDRW \
	-volname "$VOL" \
	-srcfolder "$PKG" \
	"$DMG" \
	>>$DMGLOG 2>&1

    if (( $? )); then
	rm -f "$DMG"
	echo "Error: Failed to create $DMG"
	exit 1
    fi
    
    # hdiutil attach "$DMG" -mountroot "/Volumes/$VOL" >>$DMGLOG 2>&1
    # if (( $? )); then
    # 	echo "Error: Failed to mount $DMG"
    # 	exit 1
    # fi

    # copy icons?
    
    # hdiutil detach "/Volumes/$VOL"
    # if (( $? )); then
    # 	echo "Error: Failed to unmount $DMG"
    # 	exit 1
    # fi

    # convert to compressed image, delete temp image
    echo "Compressing $VOL to $DMG..."
    hdiutil convert "$DMG" -format UDZO -o "$OUT" >>$DMGLOG 2>&1
    if (( $? )); then
	rm -f "$OUT" $DMG
	echo "Error: Failed to compress $DMG"
	exit 1
    fi

    rm -f $DMG
}

echo "Creating $NAME MacOSX disk image..."
build_dmg
echo "Created $OUT"

exit 0
