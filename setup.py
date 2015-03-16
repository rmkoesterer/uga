from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
from codecs import open
import os

execfile("uga/__version__.py")

with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'README.rst')) as f:
    long_description = f.read()

setup(
    name='uga',
	description='Universal Genome Analyst: A command line tool for analyzing genetic data',
	long_description=long_description, 
    version=version,
    url='',
    author='Ryan Koesterer',
	author_email='uga-feedback@gmail.com', 
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
    entry_points={
       'console_scripts': [
			'uga = uga.__main__:main',
           ],
       },
    packages=['uga'],
	package_data={'uga': ['data/*',]},
    classifiers = [
        'Programming Language :: Python :: 2.7',
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
		],
	zip_safe=False
)
