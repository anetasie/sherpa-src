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

import sherpa.all
import sherpa.astro.all
import sherpa.astro.ui.utils
from sherpa.utils import calc_mlr, calc_ftest, rebin, histogram1d, \
    histogram2d, gamma, lgam, erf, igamc, igam, incbet
from sherpa.data import Data1D, Data1DInt, Data2D, Data2DInt
from sherpa.astro.data import DataARF, DataRMF, DataPHA, DataIMG
from sherpa.logposterior import Prior


# We build up __all__ as we go along
__all__ = ['DataARF', 'DataRMF','DataPHA', 'DataIMG', 'Data1D', 'Data1DInt',
           'Data2D', 'Data2DInt', 'calc_mlr', 'calc_ftest', 'rebin',
           'histogram1d', 'histogram2d', 'gamma', 'lgam', 'erf', 'igamc',
           'igam', 'incbet', 'Prior']

_session = utils.Session()
_session._add_model_types(sherpa.models.basic)
_session._add_model_types(sherpa.astro.models)

if hasattr(sherpa.astro, 'xspec'):
    _session._add_model_types(sherpa.astro.xspec,
                              (sherpa.astro.xspec.XSAdditiveModel,
                               sherpa.astro.xspec.XSMultiplicativeModel))

    from sherpa.astro.xspec import get_xsabund, get_xscosmo, get_xsxsect, \
         set_xsabund, set_xscosmo, set_xsxsect
    __all__.extend(('get_xsabund', 'get_xscosmo', 'get_xsxsect',
                    'set_xsabund', 'set_xscosmo', 'set_xsxsect'))

__all__.extend(_session._export_names(globals()))


__all__ = tuple(__all__)
