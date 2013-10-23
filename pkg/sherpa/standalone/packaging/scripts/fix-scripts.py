# Copyright (c) 2009 by Enthought, Inc.
# All rights reserved.

# This script is run after core Python was installed on Unix systems.
# It does the following:
#   * update the first line in scripts in sys.prefix/bin with sys.executable
#   * update the prefix=... line in config/Makefile with sys.prefix
#
# This script is NOT supposed to be run after eggs have been installed,
# since it would be a bug in Enstaller if it does not create entry point
# scripts linking to the Python interpreter executing the Enstaller command.

import os
import sys
import re
from os.path import join, isfile, islink


bin_dir = join(sys.prefix, 'bin')
hashbang_pat = re.compile(r'#!(.+)$', re.M)


verbose = True


def fix_script(path):
    if islink(path) or not isfile(path):
        return

    fi = open(path)
    data = fi.read()
    fi.close()

    m = hashbang_pat.match(data)
    if not (m and 'python' in m.group(1).lower()):
        return

    new_data = hashbang_pat.sub('#!%s' % sys.executable, data, count=1)
    if new_data == data:
        return

    if verbose:
        print "Updating: %r" % path

    fo = open(path, 'w')
    fo.write(new_data)
    fo.close()


pre_pat = re.compile(r'^prefix=.*$', re.M)
opt_pat = re.compile(r'^OPT=.*$', re.M)
ldf_pat = re.compile(r'^LDFLAGS=.*$', re.M)
inc_pat = re.compile(r'-I\S+')
lib_pat = re.compile(r'-L\S+')
plc_pat = re.compile(r'-Wl,-R(/PLACEHOLD){5,}')

def fix_makefile(path):
    if verbose:
        print 'Fixing Mafefile: %r' % path

    fi = open(path)
    data = data_org = fi.read()
    fi.close()

    # set the current prefix
    if pre_pat.search(data):
        data = pre_pat.sub('prefix=\t\t' + sys.prefix, data)
    else:
        print "WARNING: Did not find 'prefix=...' in %r" % path

    # replace the placeholders with nothing
    data = plc_pat.sub('', data)

    # make sure the build include directories are not listed in OPT=
    m = opt_pat.search(data)
    if m:
        line = m.group()
        line = inc_pat.sub('', line).rstrip() + ' -I%s/include' % sys.prefix
        data = opt_pat.sub(line, data)
    else:
        print "WARNING: Did not find 'OPT=...' in %r" % path

    # fix LDFLAGS=
    m = ldf_pat.search(data)
    if m:
        line = m.group()
        if sys.platform == 'darwin':
            line = line.replace('-L/tmp/_py/libraries', '-L')
        else:
            line = lib_pat.sub('', line).rstrip() + ' -L%s/lib' % sys.prefix
        data = ldf_pat.sub(line, data)
    else:
        print "WARNING: Did not find 'LDFLAGS=...' in %r" % path

    # see if the data was actually changed
    if data == data_org:
        print "Hmm: Makefile is already up-to-date:"
        return

    if verbose:
        print "Updating: %r" % path

    fo = open(path, 'w')
    fo.write(data)
    fo.close()


def main():
    global verbose
    if '-q' in sys.argv or '-quiet' in sys.argv:
        verbose = False

    if verbose:
        print 'Fixing Python locations in directory: %r' % bin_dir

    for fname in os.listdir(bin_dir):
        fix_script(join(bin_dir, fname))

    from distutils.sysconfig import get_makefile_filename
    fix_makefile(get_makefile_filename())


if __name__ == '__main__':
    main()
