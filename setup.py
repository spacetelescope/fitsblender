#!/usr/bin/env python

from __future__ import division # confidence high

# We use the local copy of stsci_distutils_hack, unless
# the user asks for the stpytools version

######################################################################
# PyFITS
try:
    import pyfits
except ImportError:
    print "WARNING: PyFITS must be installed to use pywcs."
    print "         Since this is not a build-time dependency, the "
    print "         build will proceed."

import os
if os.getenv("USE_STPYTOOLS") :
    import pytools.stsci_distutils_hack as H
    pytools_version = "3.0"
else :
    import stsci_distutils_hack as H
    pytools_version = None

H.run(pytools_version = pytools_version)

