"""
Setup script for OEOMMTools

You can install oeommtools with

python setup.py install
"""

import sys,os
from os.path import relpath, join

from setuptools import setup, find_packages

if sys.argv[-1] == 'setup.py':
    print("To install, run 'python setup.py install'")
    print()

if sys.version_info[:3] < (3, 0):
    print("OEOMMTools requires Python 3.0 or later (%d.%d detected)." %
          sys.version_info[:2])
    sys.exit(-1)


descr = """
These are a collection of tools developed to integrate and mix
the OE Toolkit with the OpenMM API
"""

setup(
    name                 ='oeommtools',
    version              ='0.1.0',
    description          ='OpenEye OpenMM Tools',
    long_description     =descr,
    url                  ='https://github.com/nividic/oeommtools',
    author               ='Gaetano Calabro',
    author_email         ='gcalabro -at- eyesopen.com',
    platforms            =['Linux-64', 'Mac OSX-64', 'Unix-64'],
    packages             =find_packages()+['tests'],
    include_package_data =True,
    zip_safe             =False
)
