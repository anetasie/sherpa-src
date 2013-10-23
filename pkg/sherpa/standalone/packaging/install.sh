#!/bin/bash
echo ""
echo "Self Extracting Installer"
echo ""


# This script... do not change these
_SCRIPT_DIR=`dirname $0`
_SCRIPT_DIR=`(cd $_SCRIPT_DIR; pwd)`
_SCRIPT_NAME=`basename $0`
THIS_SCRIPT=$_SCRIPT_DIR/$_SCRIPT_NAME

# Default installation path
INSTALL_PATH=`echo $THIS_SCRIPT | sed 's_\.sh$__'`

# If BATCH_MODE=1, the script will not prompt the user and will assume defaults
BATCH_MODE=0

# Constants
PLACEHOLD=/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD/PLACEHOLD

# Prints the script usage and exits with the return code passed in
exit_with_usage()
{
    echo "usage: $0 [options]

Installs Python-PYFULLVERSION, by default to the prefix:
$INSTALL_PATH

    -h           Print this message and exit
    -b           Run in batch mode (non-interactive)
    -p <prefix>  Change the install prefix (ignoring the default)
"
    exit 2
}


# Checks the command line, exits with usage if necessary
# sets the following vars: BATCH_MODE
check_command_line()
{
    while getopts hbsp: i $*
    do
        case $i in
           h)
               exit_with_usage
               ;;
           b)
               BATCH_MODE=1
               ;;
           p)
               INSTALL_PATH="$OPTARG"
               ;;
        esac
    done
}

# Prompts the user for a directory or the default and saves the value to
# the variable $INSTALL_PATH
set_install_path()
{
    echo -n "
Python-PYFULLVERSION will be installed to this location:
$INSTALL_PATH

  * Press Enter to accept this location
  * Press CTRL-C to abort
  * or specify an alternate location.  Please ensure that your location
    contains only ASCII letters, numbers and the following punctuation
    chars: '.', '_', '-'

[$INSTALL_PATH] >>> "
    read location
    if [[ $location != "" ]]; then
        INSTALL_PATH=$location
    fi
}


# Check install path by changing to it (see if it exists) and creating the
# install path if it did not already exist.  It also checks if a file can
# be created in the directory.
check_install_path()
{
    if [ -e $INSTALL_PATH ]; then
        if [ -d $INSTALL_PATH ]; then
            echo -n "WARNING:
The install directory $INSTALL_PATH already exists.
Installing into this directory might overwrite existing files.
Are you sure you want to continue the installation? [yes|no]
[no] >>> "
            read answer
            if [[ ($answer != "yes")  && \
                  ($answer != "Yes")  && \
                  ($answer != "YES") ]]
            then
                echo "Abording installation"
                exit 2
            fi
        else
            echo "ERROR: path exists, but is not a directory: $INSTALL_PATH"
            exit 1
        fi
    fi

    mkdir -p $INSTALL_PATH
    pushd $INSTALL_PATH >/dev/null
    if (( $? )); then
        echo "ERROR: not a usable directory: $INSTALL_PATH"
        exit 1
    fi
    touch empty.txt
    if (( $? )); then
        echo "ERROR: cannot write a file into directory: $INSTALL_PATH"
        exit 1
    fi
    rm -f empty.txt
    popd >/dev/null

    # Finally, make INSTALL_PATH an absolute dir
    INSTALL_PATH=`(cd $INSTALL_PATH; pwd)`
}

# Find the start of the archive in this file and extract the archive
# in the specified install path
extract_archive()
{
    echo "Installing to $INSTALL_PATH ... please wait"

    tail -n +254 $THIS_SCRIPT | tar zxf -
    if (( $? )); then
        echo "ERROR: Could not extract tarball starting at line 245."
        exit 1
    fi

    tar xf python.tar
    rm python.tar

}

install_eggs()
{ 
    
    bin/python -E scripts/boot-enst.py eggs/enstaller*.egg
    for EGG in setuptools numpy ipython pyfits matplotlib
    do
	bin/egginst eggs/$EGG*.egg
    done
}

fix_binaries()
{

    SCRIPTS="$INSTALL_PATH/scripts"

    $SCRIPTS/chrpath -h >/dev/null
    if (( $? )); then
        echo "ERROR: Failed to run: $SCRIPTS/chrpath

One reason for this error could be that you are trying to execute a 64-bit
binary on a 32-bit system."
        exit 1
    fi

    # Fix python itself, first the Python library
    chmod u+w lib/libpythonPYVERSION.so.1.0

    $SCRIPTS/chrpath -r `pwd`/lib lib/libpythonPYVERSION.so.1.0 >/dev/null
    # now the executable(s)
    $SCRIPTS/chrpath -r `pwd`/lib bin/python  >/dev/null
    $SCRIPTS/chrpath -r `pwd`/lib bin/pythonPYVERSION  >/dev/null

    # Now we have Python working, make sure we can call interpreter
    echo "Trying to run Python interpreter:"
    bin/python -E -c "import sys; print sys.prefix"
    if (( $? )); then
        echo "ERROR: Cannot run Python interpreter: $INSTALL_PATH/bin/python"
        exit 1
    fi

}

fix_rpaths()
{
    
    echo "Fixing RPATH in binaries..."
    chmod -R u+w lib

    $INSTALL_PATH/bin/python -E scripts/fix-rpath.py
    $INSTALL_PATH/bin/python -E scripts/fix-scripts.py -q
    echo "Compiling all Python modules..."
    $INSTALL_PATH/bin/python -E $INSTALL_PATH/lib/pythonPYVERSION/compileall.py -f -q \
        -x 'bad_coding|badsyntax|py3_|custom_tools' $INSTALL_PATH/lib/pythonPYVERSION

}

cleanup()
{

    /bin/rm -rf scripts eggs EGG-INFO

}


# Print a message indicating where the new install is
print_goodbye()
{
    echo "done.

    As the last step, you should edit your .bashrc or prepend
    the Python-PYFULLVERSION install path:

        $INSTALL_PATH/bin

    Thank you for installing Python-PYFULLVERSION!
"
    if [[ $PYTHONPATH != "" ]] 
    then
        echo "
WARNING: the PYTHONPATH environment variable is currently set.  This might
         cause unexpected behavior when running the newly installed Python
         interpreter.
"
    fi
}

# INSTALL SCRIPT
check_command_line $*

if (( ! $BATCH_MODE ))
then
    set_install_path
fi

check_install_path

pushd $INSTALL_PATH >/dev/null
extract_archive
fix_binaries
install_eggs
fix_rpaths
cleanup

popd >/dev/null

print_goodbye

# End of script, everything after this must be ignored by the shell
# so we call exit explicitly
exit 0

# TARFILE CONCATENATED BELOW TAG (tag includes a single \n at the end)
__ARCHIVE_BELOW__
