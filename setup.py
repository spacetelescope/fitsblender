#!/usr/bin/env python

from setuptools import find_packages
from setuptools import setup
from pathlib import Path


this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

PACKAGENAME = 'fitsblender'

setup(
    name=PACKAGENAME,
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    install_requires=[
        'astropy>=5.0.4',
        'numpy',
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
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/spacetelescope/fitsblender',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
