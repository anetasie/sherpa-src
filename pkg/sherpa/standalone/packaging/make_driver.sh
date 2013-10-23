#!/bin/bash

SCRIPT_DIR=`dirname $0`
SCRIPT_DIR=`(cd $SCRIPT_DIR; pwd)`
SCRIPT_NAME=`basename $0`
THIS_SCRIPT=$SCRIPT_DIR/$SCRIPT_NAME
VERSION_STRING="$SCRIPT_NAME v1.0"

PYTHONVER="2.7.2"

usage()
{
    echo "${VERSION_STRING}"
    echo 
    echo "Usage: ${SCRIPT_NAME} [options...]"
    echo
    echo "  -h --help               Print this message"
    echo "  -v --version            Python version to build e.g 2.7.2"
    echo "  -b --build              System to build { Linux64 Linux OSXI OSX64 }"
    echo "  --without-installer     Build Python but omit the installer"

}

# Determine the default build type

ARCH=`arch`
KERNEL=`uname -s`
BUILD=
WITHOUT=

case "$KERNEL" in

    Linux* )
	ARCH=`uname -m`
	if [ "$ARCH" != "x86_64" ] ; then
	    BUILD="Linux"
	else
	    BUILD="Linux64"
	fi
	;;

    Darwin* | darwin* )
	ARCH=`uname -p`
	case "$ARCH" in
	    i386* )
		if [ "x`sysctl hw.cpu64bit_capable | grep 1`" != x ] ; then
		    BUILD="OSX64"
		else
		    BUILD="OSXI"
		fi
		;;
	    
	    x86_64* )
		BUILD="OSX64"
		;;

	    * )
		echo "ERROR: $KERNEL $ARCH not supported" >&2
		exit 1
		;;
	esac
	;;

    * )
	echo "ERROR: $KERNEL not supported" >&2
	exit 1
	;;
esac

SYSTEM="$BUILD"


# Utilize command line switches a la configure

ac_prev=
ac_dashdash=
for ac_option in $@
do
    # If the previous option needs an argument, assign it.
    if test -n "$ac_prev"; then
	eval $ac_prev=\$ac_option
	ac_prev=
	continue
    fi

    case $ac_option in
	*=*)	ac_optarg=`expr "X$ac_option" : '[^=]*=\(.*\)'` ;;
	*)	ac_optarg=yes ;;
    esac
    
    case $ac_dashdash$ac_option in
	--)
	    ac_dashdash=yes ;;
	
	-build | --build | -b )
	    ac_prev=BUILD ;;
	
	-build=* | --build=* | --buil=* | --bui=* | --bu=* ) 
            BUILD=$ac_optarg ;;

	-help | --help | -h* | --h* )
	    usage
	    exit 0 ;;

	--version | -v )
	    ac_prev=PYTHONVER ;;

	-version | --version | -v* | --v* )
	    PYTHONVER=$ac_optarg ;;
	
	--without-i* | -without-i* )
	    WITHOUT=yes
	    ;;

	-* | --* | * )
	    echo "Error: unrecognized option $ac_option" >&2
	    exit 1 ;;
	
    esac
done


# Check the user's choice against the current architecture
case $BUILD in

    OSX* )
	if [ $KERNEL = "Linux" ]; then
	    echo "ERROR: Trying to build $BUILD on $KERNEL"
	    exit 1
	fi
	;;

    Linux* )
	if [ $KERNEL = "Darwin" ]; then
	    echo "ERROR: Trying to build $BUILD on $KERNEL"
	    exit 1
	fi
	;;

    * )
	echo "ERROR: unknown build platform $BUILD"
	exit 1
	;;
esac


VERSION=$PYTHONVER
MAJOR=`echo $VERSION | awk '{print substr($1,1,3)}'`
PYTHON="Python-$VERSION"
SRC=$SCRIPT_DIR
BUILDPY="$SRC/make_python.sh"
BUILDOTS="$SRC/make_ots.sh"
INSTALLSCRIPT="$SRC/install.sh"
INSTALLER="$SRC/make_installer.sh"
PACKAGEBUILDER="$SRC/make_python_packages.sh"

# Build Python OTS

$BUILDOTS $BUILD $PYTHON
if (( $? )); then
    echo "Python OTS build failed!"
    exit 1
fi

# Install distribution helper scripts

SCRIPTDIR="$SRC/$BUILD/$PYTHON/INSTALL/scripts"
mkdir -p $SCRIPTDIR 
cp $SRC/scripts/* $SCRIPTDIR

case $BUILD in

    Linux* )
	CHRPATH="$SRC/$BUILD/$PYTHON/INSTALL/bin/chrpath"
	if [ -f $CHRPATH ]; then
	    cp $CHRPATH $SCRIPTDIR
	fi
	;;

    OSX* )
	INSTALL_NAME_TOOL="/usr/bin/install_name_tool"
	if [ -f $INSTALL_NAME_TOOL ]; then
	    cp -f $INSTALL_NAME_TOOL $SCRIPTDIR
	fi

	INSTALLSCRIPT="$SRC/install.sh.mac"
	;;
esac


# if [ ! -f $SRC/$BUILD/$PYTHON/INSTALL/bin/python ]; then
#     # Clean out the OTS executables in bin/
#     /bin/rm -f $SRC/$BUILD/$PYTHON/INSTALL/bin/*
# fi


# Build Python

$BUILDPY $BUILD $VERSION
if (( $? )); then
    echo "Python build failed!"
    exit 1
fi

SYS=`uname -s`

# Sanitize the dylibs on Mac
if [ "$SYS" = "Darwin" ]; then

    cd $SRC/$BUILD/$PYTHON/INSTALL/lib
    DYLIBS=`find . -print -maxdepth 1 | grep dylib`
    for ii in $DYLIBS
    do
	echo "Fixing dylib: $ii"
	chmod +w $ii
	NAME=`echo $ii | awk '{print substr($1,3)}'`
	install_name_tool -id "@rpath/$NAME" $ii
    done
    cd $SRC
fi

# Save clean version of Python off to the side for installer
TARFILE="python.tar"
TAR=tar
if [ "$SYS" = "Darwin" ]; then
    TAR=gnutar
fi

tar_python()
{
    DIR="$BUILD/$PYTHON/INSTALL"
    cd $DIR
    if [ ! -f ../$TARFILE ]; then
	echo "Tarring $PYTHON ..."
	$TAR -cf ../$TARFILE bin docs include lib scripts share
	if (( $? )); then
	    echo "Error: Tarring $DIR python failed!"
	    exit 1
	fi
    else
	echo "$PYTHON tar file found, skipping archive..."
    fi
    cd $SRC
}
    
tar_python


# make python packages
$PACKAGEBUILDER $BUILD $MAJOR $VERSION
if (( $? )); then
    echo "Build failed!"
    exit 1
fi

# make installer
if [ "x$WITHOUT" = "x" ]; then
    cd $BUILD/$PYTHON
    $INSTALLER INSTALL "$PYTHON-$BUILD" $MAJOR $VERSION $INSTALLSCRIPT
    if (( $? )); then
	echo "Installer failed!"
	exit 1
    fi
fi
