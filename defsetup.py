#!/usr/bin/env python

## this is obsolete - see setup.cfg

from __future__ import division # confidence high

CONTACT = "Michael Droettboom"
EMAIL = "support@stsci.edu"
VERSION = '0.1'

from distutils.core import setup, Extension
from os.path import join
import os.path
import sys


######################################################################
# DISTUTILS SETUP
libraries = []

pkg = ["fitsblender", "fitsblender.tests"]

setupargs = {
    'version' :     VERSION,
    'description':  "Aggregate values in FITS headers",
    'author' :      CONTACT,
    'author_email': EMAIL,
    'url' :         "",
    'platforms' :   ["unix","windows"],
    'package_dir' : {pkg[0]: 'lib/fitsblender', pkg[1]: 'lib/tests'},
    'data_files' : [
        ( 'fitsblender/pars',['lib/fitsblender/pars/*']),
        ( 'fitsblender/tests', ['lib/tests/*.fits']),
        ( 'fitsblender',['lib/fitsblender/*.rules']),
        ],
}

