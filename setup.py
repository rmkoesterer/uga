from __future__ import print_function
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import codecs
import os
import sys
import re

here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    # intentionally *not* adding an encoding option to open
    return codecs.open(os.path.join(here, *parts), 'r').read()

VERSIONFILE="uga/version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))	

long_description = read('README.rst')

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['--strict', '--verbose', '--tb=long', 'tests']
        self.test_suite = True

    def run_tests(self):
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)

setup(
    name='uga',
    version=verstr,
    url='',
    license='GPL2.0',
    author='Ryan Koesterer',
    tests_require=['pytest'],
    install_requires=['biopython==1.64', 
						'multi-key-dict==2.0.1', 
						'numpy==1.9.1', 
						'pandas==0.15.1', 
						'patsy==0.3.0', 
						'progressbar==2.3', 
						'psutil==2.1.3', 
						'pytabix==0.1', 
						'python-dateutil==2.2', 
						'pytz==2014.9', 
						'rpy2==2.5.2', 
						'scipy==0.14.0', 
						'six==1.8.0', 
						'statsmodels==0.6.0', 
						'virtualenv==1.10.1', 
						'wsgiref==0.1.2',
						'singledispatch==3.4.0.3'],
    cmdclass={'test': PyTest},
    author_email='rmkoesterer@gmail.com',
    description='Universal Genome Analyst',
    long_description=long_description,
    entry_points={
       'console_scripts': [
           'uga = uga.__main__:main',
           ],
       },
	#scripts=['uga/'],
    packages=['uga'],
    #include_package_data=True,
    platforms='any',
    #test_suite='sandman.test.test_sandman',
    #package_data={'uga': ['uga/*']},
    classifiers = [
        'Programming Language :: Python :: 2',
        'Development Status :: Pre-release',
        'Natural Language :: English',
        'Environment :: Genetics Data Analysis',
        'Intended Audience :: Users',
        'License :: GPL2.0',
        'Operating System :: Unix',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Software Development :: Libraries :: Application Frameworks',
        ],
    extras_require={
        'testing': ['pytest'],
      }
)