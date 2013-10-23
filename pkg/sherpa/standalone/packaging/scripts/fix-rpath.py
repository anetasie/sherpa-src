# Copyright (c) 2009 by Enthought, Inc.
# All rights reserved.

"""
Update the RPATH settings in the ELF headers of various executable binaries
and libraries installed within this Python environment.
"""

import os
import sys
from os.path import join, islink, isfile


placeholder = '/PLACEHOLD' * 20


def is_elf(filepath):
    """
    Return True if the specified file is an ELF binary file.
    """
    if filepath.endswith(('.py', '.pyc', '.pyo', '.h', '.a')):
        return False

    if islink(filepath) or not isfile(filepath):
        return False

    fi = open(filepath, 'rb')
    head = fi.read(4)
    fi.close()

    return head == '\x7fELF'


def set_rpath(filepath, rpath):
    """
    Replace the placeholder with the given rpath in filepath
    """
    def replace(data, old, new):
        """
        Replace the bytes `old` with `new` and return the new data
        """
        assert len(new) == len(old)
        return data.replace(old, new)

    fi = open(filepath, 'rb')
    data = fi.read()
    fi.close()

    # First try to replace the NULL terminated placeholder with the new rpath
    # padded with zeros
    new_data = replace(
        data,
        placeholder + chr(0),
        rpath + (len(placeholder) - len(rpath) + 1) * chr(0))

    # Now try to replace the placeholder itself.  This will be used when
    # the placeholder is followed by another path.
    new_data = replace(
        new_data,
        placeholder,
        (rpath + ":/" + len(placeholder) * 'A')[:len(placeholder)])

    if new_data != data:
        fo = open(filepath, 'wb')
        fo.write(new_data)
        fo.close()
        #print 'Fixed RPATH in %s' % filepath
        sys.stdout.write('.')
        sys.stdout.flush()


def set_rpath_r(path):
    """
    Update the RPATH in all ELF files within a directory tree.
    """
    targets = [join(sys.prefix, 'lib')]

    # Fatal error if we can't set the desired paths into the space we have
    # available!
    rpaths = os.pathsep.join(targets)
    if len(rpaths) > len(placeholder):
        raise Exception("ERROR: New rpath exceeds %i characters" %
                        len(placeholder))

    # Walk the specified directory, updating all the binaries along the way.
    for root, dirs, files in os.walk(path):
        for fn in files:
            fpath = join(root, fn)
            if is_elf(fpath):
                set_rpath(fpath, rpaths)


def main():
    """
    Fix up all ELF binaries in the bin and lib directory.
    """
    print 'Updating RPATH in ELF headers'
    for subdir in ['bin', 'lib']:
        path = join(sys.prefix, subdir)
        #print 'Updating RPATH in ELF headers in %s' % path
        set_rpath_r(path)
    print


if __name__ == '__main__':
    main()
