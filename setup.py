#!/usr/bin/env python
import relic.release
from setuptools import setup, find_packages


version = relic.release.get_info()
relic.release.write_template(version, 'lib/fitsblender')


setup(
    name = 'fitsblender',
    version = version.pep386,
    author = 'Michael Droettboom',
    author_email = 'help@stsci.edu',
    description = 'Aggregate values in FITS headers',
    url = 'https://github.com/spacetelescope/fitsblender',
    classifiers = [
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    install_requires = [
        'astropy',
        'nose',
        'numpy',
        'pytools',
        'sphinx',
        'stsci.sphinxext',
        'stsci.tools',
    ],

    package_dir = {
        '':'lib'
    },
    packages = find_packages('lib'),
    package_data = {
        'fitsblender': [
            'pars/*',
            'tests/*.fits',
            '*.rules',
            'LICENSE.txt'
        ]
    },
)
