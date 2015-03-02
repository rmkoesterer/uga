from __future__ import print_function
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import codecs
import os
import sys
import re

here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    return codecs.open(os.path.join(here, *parts), 'r').read()

VERSIONFILE="uga/__version__.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))	

long_description = read('README.rst')

setup(
    name='uga',
    version=verstr,
    url='',
    license='The MIT License: http://www.opensource.org/licenses/mit-license.php',
    author='Ryan Koesterer',
	install_requires=['singledispatch', 
						'rpy2', 
						'multi-key-dict', 
						'numpy', 
						'pandas', 
						'progressbar', 
						'psutil', 
						'pytabix', 
						'scipy', 
						'biopython', 
						'plinkio', 
						'pysam', 
						'PyVCF', 
						'plinkio'], 
    author_email='rmkoesterer@gmail.com',
    description='Universal Genome Analyst',
    long_description=long_description,
    entry_points={
       'console_scripts': [
           'uga = uga.__main__:main',
		   'quga = uga.QUga:main',
           ],
       },
    packages=['uga'],
    platforms='any',
    classifiers = [
        'Programming Language :: Python :: 2',
        'Development Status :: Pre-release',
        'Natural Language :: English',
        'Environment :: Genetics Data Analysis',
        'Intended Audience :: Users',
        'License :: The MIT License',
        'Operating System :: Unix',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Software Development :: Libraries :: Application Frameworks',
        ]
)