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
	install_requires=['singledispatch==3.4.0.3', 
						'rpy2==2.5.2', 
						'multi-key-dict==2.0.1', 
						'numpy==1.9.1', 
						'pandas==0.15.1', 
						'progressbar==2.3', 
						'psutil==2.1.3', 
						'pytabix==0.1', 
						'scipy==0.14.0', 
						'biopython==1.64'], 
    author_email='rmkoesterer@gmail.com',
    description='Universal Genome Analyst',
    long_description=long_description,
    entry_points={
       'console_scripts': [
           'uga = uga.__main__:main',
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