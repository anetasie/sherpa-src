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

import numpy
from sherpa.utils.err import InstrumentErr
from sherpa.models.model import ArithmeticFunctionModel, NestedModel
from sherpa.astro.models import MultiResponseSumModel

__all__ = ('standard_fold', 'multiresponse_fold', 'pileup_fold')


def standard_fold(data, model):
    if data.exposure:
        model = data.exposure * model

    # clear out any previous response filter
    data.notice_response(False)

    arf, rmf = data.get_response()

    # Notice response energy bins according to noticed channels
    # Noticing at fold time ensures that the mask is finalized before
    # fitting
    if numpy.iterable(data.mask):
        data.notice_response(True)

    if arf is not None:
        args = ()

        # PHA <=> ARF case
        if data.bin_lo is not None and data.bin_hi is not None:
            if len(data.bin_lo) != len(arf.energ_lo):
                args = data._get_ebins(group=False)

        # RMF <=> ARF case
        if rmf is not None:
            if len(arf.energ_lo) != len(rmf.energ_lo):
                args = rmf.get_indep()

        model = model.apply(arf.apply_arf, arf.get_indep(), args)

    if rmf is not None:
        args = ()

        # PHA <=> RMF case
        if arf is None and data.bin_lo is not None and data.bin_hi is not None:
            if len(data.bin_lo) != len(rmf.energ_lo):
                args = data._get_ebins(group=False)

        model = model.apply(rmf.apply_rmf, rmf.get_indep(), args)
    return model


def multiresponse_fold(data, model):
    models = []

    # clear out any previous response filter
    data.notice_response(False)

    # clear out the previous hi-res ARF grid
    data._multi_resp_grid=None

    # Notice response energy bins according to noticed channels
    # Noticing at fold time ensures that the mask is finalized before
    # fitting
    if numpy.iterable(data.mask):
        data.notice_response(True)

    #for arf, rmf in data._responses.values():
    for id in data.response_ids:
        arf, rmf = data.get_response(id)
        m = None
        if arf is not None:
            m = ArithmeticFunctionModel(arf.apply_arf)
        if rmf is not None:
            if m is not None:
                m = NestedModel( rmf.apply_rmf,  arf.apply_arf)
            else:
                m = ArithmeticFunctionModel(rmf.apply_rmf)
        models.append(m)

    model = MultiResponseSumModel(models, model)

    if data.exposure:
        model = data.exposure * model

    return model


def pileup_fold(data, model, pileupmodel):
    # clear out any previous response filter
    data.notice_response(False)

    arf, rmf = data.get_response()
    err_msg = None
    
    if arf is None:
        err_msg = 'does not have an associated ARF'
    elif data.exposure is None:
        err_msg = 'does not specify an exposure time'

    if err_msg:
        raise InstrumentErr('baddata', data.name, err_msg)

    # Currently, the response is NOT noticed using pileup

    # ARF convolution done inside ISIS pileup module
    # on finite grid scale
    #model = model.apply(arf.apply_arf)
  
    model = model.apply(pileupmodel, data.exposure, arf.energ_lo, arf.energ_hi,
                        arf.specresp, model)
    if rmf is not None:
        model = model.apply(rmf.apply_rmf)

    return model
