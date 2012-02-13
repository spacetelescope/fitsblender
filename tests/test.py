# Copyright (C) 2008-2010 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

import glob
import os

import numpy as np
import pyfits
from numpy.testing import \
    assert_equal, assert_almost_equal, assert_array_almost_equal, assert_raises

import fitsblender

ROOT_DIR = None
def setup():
    global ROOT_DIR
    ROOT_DIR = os.path.dirname(__file__)

def _test_fitsblender(files):
    spec = [
        ('CRPIX1', 'crpix', np.average),
        ('crpix1', 'crpix_lower', np.average),
        ('CRPIX1', 'crpix_constant', np.average, 'constant', 0),
        ('CRPIX1', 'crpix'),
        ('CRPIX1', 'crpix_all', np.array),
        ('CRPIX1', 'crpix_range', (np.min, np.max)),
        ('OBSERVER', 'observer'),
        ('CRPIX1', 'crpix_sum', np.sum),
        ('CRPIX1', 'crpix_last', lambda x: x[-1]),
        ('CRPIX1', 'crpix_stddev', np.std),
        ]

    d, t = fitsblender.fitsblender(files, spec)

    print d
    print t
    print t.dtype

    assert_equal(d['crpix'], d['crpix_lower'])
    assert_almost_equal(d['crpix'], 77.6666666666666)
    assert_almost_equal(d['crpix_constant'], 51.7777777777777777)
    assert_array_almost_equal(d['crpix_all'], [1, 1, 1, 256.5, -213.5, 420.])
    assert_array_almost_equal(t['crpix'], [np.nan, 1, 1, 1, np.nan, np.nan, 256.5, -213.5, 420.])
    assert_array_almost_equal(d['crpix_range'], [-213.5, 420.0])
    assert t['observer'][0] == 'A. A. Zdziarski'
    assert_almost_equal(d['crpix_sum'], 466.0)
    assert_almost_equal(d['crpix_last'], 420.0)
    assert_almost_equal(d['crpix_stddev'], 204.77012857239592)

def test_filenames():
    files = [
        (x, 0) for x in glob.glob(os.path.join(ROOT_DIR, "*.fits"))]
    files.sort()

    _test_fitsblender(files)

def test_header_objects():
    files = glob.glob(os.path.join(ROOT_DIR, "*.fits"))
    files.sort()
    headers = [pyfits.open(x)[0].header for x in files]

    _test_fitsblender(headers)

def test_nospec():
    """
    Test that no header list and no spec doesn't raise an exception.
    """
    d, t = fitsblender.fitsblender([], [])
    assert d == {}
    assert len(t) == 0

def test_raise():
    def raises():
        d, t = fitsblender.fitsblender(files, spec)

    spec = [
        ('CRPIX1', 'crpix_raise', np.average, 'raises')
        ]

    files = [
        (x, 0) for x in glob.glob(os.path.join(ROOT_DIR, "*.fits"))]
    files.sort()

    assert_raises(ValueError, raises)

def test_raise_table():
    def raises():
        d, t = fitsblender.fitsblender(files, spec)

    spec = [
        ('CRPIX1', 'crpix_raise', None, 'raises')
        ]

    files = [
        (x, 0) for x in glob.glob(os.path.join(ROOT_DIR, "*.fits"))]
    files.sort()

    assert_raises(ValueError, raises)
