#!/usr/bin/env python

from setuptools import find_packages
from setuptools import setup

PACKAGENAME = 'fitsblender'

setup(
    name=PACKAGENAME,
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    install_requires=[
        'astropy',
        'numpy',
        'six',
        'stsci.tools',
    ],
    extras_require={
        'docs': [
            'sphinx',
            'numpydoc',
        ],
        'test': [
            'pytest',
            'pytest-cov',
            'codecov',
        ],
    },
    packages=find_packages(),
    package_data={
        PACKAGENAME: [
            '*.rules',
            'pars/*',
            'tests/*.fits',
        ],
    },
    author='Michael Droettboom',
    author_email='help@stsci.edu',
    description='Aggregate values in FITS headers',
    url='https://github.com/spacetelescope/fitsblender',
    classifiers = [
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
