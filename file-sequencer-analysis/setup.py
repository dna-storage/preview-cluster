# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.core import setup, Extension

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()
   
setup(name='generate', version='1.0', ext_modules=[Extension('generate',
                                                            ['random_int.c'])])

setup(
    name='dnastorage',
    version='0.1.0',
    description='DNA-based data storage modeling and simulation package',
    long_description=readme,
    author='James Tuck',
    author_email='jtuck@ncsu.edu',
    url='https://github.ncsu.edu/jtuck/',
    license=license,
    packages=find_packages(exclude=( 'tests', 'docs', 'tools',
                                    'other_software'))
)
