#!/usr/bin/env python
# 
#  Copyright (C) 2011  Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#


import os
import sys
import stat
import shutil


perm = stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH

def permissions_backup(filename, bakname):
    os.chmod(filename, perm)
    if os.path.isfile(bakname):
        os.chmod(bakname, perm)
        shutil.copy(bakname, filename)
    else:
        shutil.copy(filename, bakname)


def write(fd, line):

    fd.seek(0)
    fd.write(line)
    fd.flush()
    fd.close()


def update_xspec(root):

    initfile = os.path.join(root, 'sherpa', 'astro', 'xspec', '__init__.py')
    xspecfile = os.path.join(root, 'sherpa', 'astro', 'xspec', 'src', '_xspec.cc')

    permissions_backup(initfile, initfile+'.bak')
    permissions_backup(xspecfile, xspecfile+'.bak')

    fd = open(initfile, "rw+")
    line = fd.read()

    xspec__init__ = """
import os
from sherpa.astro.xspec._xspec import xsinit
if os.environ.get('HEADAS', None) is None:
    xsinit(os.path.join(os.path.dirname(__file__), 'spectral'))
else:
    xsinit(os.environ.get('HEADAS'))
import string
    """
    line = line.replace("import string", xspec__init__)
    write(fd, line)


    fd = open(xspecfile, "rw+")
    line = fd.read()

    line = line.replace("void FNINIT(void);", "void FNINIT(char* headas);")
    line = line.replace("int _sherpa_init_xspec_library();", "int _sherpa_init_xspec_library(char* headas=0x0);")
    line = line.replace("int _sherpa_init_xspec_library()", "int _sherpa_init_xspec_library(char* headas)")
    line = line.replace("FNINIT();", "FNINIT(headas);")

    xspec_xspec_func = """
static PyObject* _init( PyObject *self, PyObject *args )
 { 
   char* headas = NULL;

   if ( !PyArg_ParseTuple( args, (char*)"s", &headas ) )
      return NULL;

   if ( EXIT_SUCCESS != _sherpa_init_xspec_library(headas) )
      return NULL;

   Py_RETURN_NONE;

 }

static PyMethodDef XSpecMethods[] = {
 { (char*)"xsinit", (PyCFunction)_init, METH_VARARGS,
   (char*)"XSPEC initialization, HEADAS." },
"""

    xspec_xspec_decl_target = "static PyMethodDef XSpecMethods[] = {"

    line = line.replace(xspec_xspec_decl_target, xspec_xspec_func)
    write(fd, line)


def update_rcfile(root, oldversion, newversion):

    rcfile = os.path.join(root, 'sherpa', 'sherpa.rc')
    permissions_backup(rcfile, rcfile+'.bak')

    fd = open(rcfile, "rw+")

    line = fd.read()

    line = line.replace(oldversion, newversion)
    line = line.replace(": chips", ": pylab")
    line = line.replace(": crates", ": pyfits")

    write(fd, line)


def update_version(root, oldversion, newversion):

    setupfile = os.path.join(root, 'setup.py')
    initfile = os.path.join(root, 'sherpa', '__init__.py')

    permissions_backup(setupfile, setupfile+'.bak')
    permissions_backup(initfile, initfile+'.bak')

    fd = open(setupfile, "rw+")
    line = fd.read()
    line = line.replace(oldversion, newversion)

    write(fd, line)

    fd = open(initfile, "rw+")
    line = fd.read()
    line = line.replace(oldversion, newversion)

    oldint = oldversion.replace('.', '0')
    newint = newversion.replace('.', '0')

    line = line.replace(oldint, newint)

    write(fd, line)
   


def main():
    if len(sys.argv) < 2:
        raise RuntimeError("Please provide a root argument and old and new version arguments")

    root = str(sys.argv[1]).strip()
    old = str(sys.argv[2]).strip()
    new = str(sys.argv[3]).strip()

    update_xspec(root)
    update_rcfile(root, old, new)
    update_version(root, old, new)

if __name__ == "__main__":
    main()
