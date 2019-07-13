from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
from codecs import open
import os
from Cython.Build import cythonize
import numpy as np
from uga import __version__
import subprocess

with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'README.rst')) as f:
    long_description = f.read()

githash = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip()

if githash != "":
	__version__ = __version__ + '-' + githash

with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'uga/__version__.py'), "w") as f:
	f.write('version="' + __version__ + '"')

setup(
    name='uga',
	description='Universal Genome Analyst (uga) is a tool designed to assist biomedical researchers in complex genomic data analysis',
	long_description=long_description, 
    version=__version__,
    url='',
    author='Ryan Koesterer',
	author_email='uga-feedback@gmail.com', 
	ext_modules = cythonize(["uga/Geno.pyx","uga/Geno.pxd","uga/Model.pyx","uga/Variant.pyx","uga/Variant.pxd"]), 
	install_requires=['singledispatch', 
						'rpy2', 
						'numpy', 
						'pandas', 
						'progressbar', 
						'psutil', 
						'biopython', 
						'pysam', 
						'Cython', 
						'scipy'], 
    entry_points={
       'console_scripts': [
			'uga = uga.__main__:main',
           ],
       },
    packages=['uga'], 
	include_dirs = [np.get_include()], 
	package_data={'uga': ['settings.ini',]},
    classifiers = [
        'Programming Language :: Python :: 2.7',
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
		],
	zip_safe=False
)
