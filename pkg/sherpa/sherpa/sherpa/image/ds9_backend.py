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

from itertools import izip
import numpy
import time
import DS9
from os import access, R_OK
from sherpa.utils import get_keyword_defaults
from sherpa.utils.err import DS9Err

imager = DS9.DS9Win(DS9._DefTemplate, False)


def close():
    if imager.isOpen():
        imager.xpaset("quit")

def delete_frames():
    if not imager.isOpen():
        raise DS9Err('open')
    try:
        imager.xpaset("frame delete all")
        return imager.xpaset("frame new")
    except:
        raise DS9Err('delframe')
        
def get_region():
    if not imager.isOpen():
        raise DS9Err('open')
    try:
        imager.xpaset("regions format ciao")
        imager.xpaset("regions strip yes")
        return imager.xpaget("regions")
    except:
        raise DS9Err('retreg')
    
def image(arr, newframe=False, tile=False):
    if not imager.isOpen():
        imager.doOpen()
    if newframe is True:
        try:
            imager.xpaset("frame new")
            imager.xpaset("frame last")
        except:
            raise DS9Err('newframe')
    try:
        if tile is True:
            imager.xpaset("tile yes")
        else:
            imager.xpaset("tile no")
    except:
        raise DS9Err('settile')
    time.sleep(1)
    try:
        imager.showArray(arr)
    except:
        raise DS9Err('noimage')

def _set_wcs(keys):
    eqpos, sky, name = keys

    phys = ''
    wcs  = "OBJECT = '%s'\n" % name

    if eqpos is not None:
        wcrpix  = eqpos.crpix
        wcrval  = eqpos.crval
        wcdelt  = eqpos.cdelt

    if sky is not None:
        pcrpix = sky.crpix
        pcrval = sky.crval
        pcdelt = sky.cdelt

        # join together all strings with a '\n' between each
        phys = '\n'.join(["WCSNAMEP = 'PHYSICAL'",
                          "CTYPE1P = 'x       '",
                          'CRVAL1P = %.14E' % pcrval[0],
                          'CRPIX1P = %.14E' % pcrpix[0],
                          'CDELT1P = %.14E' % pcdelt[0],
                          "CTYPE2P = 'y       '",
                          'CRVAL2P = %.14E' % pcrval[1],
                          'CRPIX2P = %.14E' % pcrpix[1],
                          'CDELT2P = %.14E' % pcdelt[1]])

        if eqpos is not None:
            wcdelt = wcdelt * pcdelt
            wcrpix = ((wcrpix - pcrval) /
                      pcdelt + pcrpix )

    if eqpos is not None:
        # join together all strings with a '\n' between each
        wcs = wcs + '\n'.join(["RADECSYS = 'ICRS    '",
                               "CTYPE1  = 'RA---TAN'",
                               'CRVAL1  = %.14E' % wcrval[0],
                               'CRPIX1  = %.14E' % wcrpix[0],
                               'CDELT1  = %.14E' % wcdelt[0],
                               "CTYPE2  = 'DEC--TAN'",
                               'CRVAL2  = %.14E' % wcrval[1],
                               'CRPIX2  = %.14E' % wcrpix[1],
                               'CDELT2  = %.14E' % wcdelt[1]])

    # join the wcs and physical with '\n' between them and at the end
    return ('\n'.join([wcs,phys]) + '\n')

def wcs(keys):

    if not imager.isOpen():
        raise DS9Err('open')

    info = _set_wcs( keys )

    try:
        # use stdin to pass the WCS info
        imager.xpaset('wcs replace', info)
    except:
        raise DS9Err('setwcs')


def open():
    imager.doOpen()

def set_region(reg):
    if not imager.isOpen():
        raise DS9Err('open')
    try:
        imager.xpaset("regions format ciao")
    except:
        raise DS9Err('setreg')
    try:
        if (access(reg, R_OK) is True):
            imager.xpaset("regions load " + "'" + reg + "'")
        else:
            imager.xpaset("regions command " + "'" + reg + "'")
    except:
        raise DS9Err('badreg', str(reg))

def xpaget(arg):
    if not imager.isOpen():
        raise DS9Err('open')
    return imager.xpaget(arg)

def xpaset(arg):
    if not imager.isOpen():
        raise DS9Err('open')
    return imager.xpaset(arg)
