# 
#  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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

"""

Modeling and fitting package for scientific data analysis

Sherpa is a modeling and fitting package for scientific data analysis.
It includes general tools of interest to all users as well as
specialized components for particular disciplines (e.g. astronomy).

Note that the top level sherpa package does not import any
subpackages.  This allows the user to import a particular component
(e.g. sherpa.optmethods) without pulling in any others.  To import all
the standard subpackages, use 'import sherpa.all' or
'from sherpa.all import *'.

"""


import logging
import os
import os.path
import sys
from sherpa.utils import SherpaTest


__all__ = ('banner', 'get_include', 'test')

__versionstr__ = '4.2.1'
__version__ = 40201

class Formatter(logging.Formatter):
    def format(self, record):
        if record.levelno > logging.INFO:
            msg = '%s: %s' % (record.levelname, record.msg)
        else:
            msg = record.msg
        return msg

log = logging.getLogger('sherpa')
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(Formatter())
log.addHandler(handler)
log.setLevel(logging.INFO)

del Formatter, log, handler


_banner = """
-----------------------------------------------------
Welcome to Sherpa: CXC's Modeling and Fitting Package
-----------------------------------------------------
Version: CIAO %s
""" % __versionstr__


def banner(file=sys.stdout):
    "Print a welcome message to the specified file object"
    print >> file, _banner


def get_include():
    "Get the root path for installed Sherpa header files"
    
    return os.path.join(os.path.dirname(__file__), 'include')

def get_config():
    "Get the path for the installed Sherpa configuration file"

    home_dir = None
    config = None

    # If NOSHERPARC is set, read in system config file
    # ignore any user config file
    if (os.environ.has_key('NOSHERPARC') == True):
        return os.path.join(os.path.dirname(__file__), 'sherpa.rc')
    
    # If SHERPARC is set, read in config file from there,
    # and ignore default location
    if (os.environ.has_key('SHERPARC') == True):
        config = os.environ.get('SHERPARC')
        if os.path.isfile(config):
            return config
        
    # SHERPARC was not set, so look for .sherpa.rc in default
    # location, which is user's home directory.
    home_dir = os.environ.get('HOME')
    config = os.path.join(home_dir, '.sherpa.rc')

    if os.path.isfile(config):
        return config

    # If no user config file is set, fall back to system config file
    return os.path.join(os.path.dirname(__file__), 'sherpa.rc')

    
def test(level=1, verbosity=1, datadir=None):
    """
    
    Run the Sherpa test suite, testing all available subpackages
    (including the discipline-specific ones)
    
    """

    import sherpa.all
    import sherpa.astro.all
    SherpaTest().test(level, verbosity, datadir)
