Universal Genome Analyst
************************

Universal Genome Analyst (**uga**) is an open, flexible, and efficient tool for the distribution, management, and visualization of whole genome data analyses. 
It is designed to assist biomedical researchers in complex genomic data analysis through the use of a low level interface between the powerful R statistical environment and Python, allowing
for rapid integration of emerging analytical strategies. This project uses `Cython`_ for a significant reduction in computation time and researchers with access to a high performance computing cluster or 
with access to multiple cores will find time-saving features for parallel analysis using a flexible, yet controlled, commandline interface.

This software is currently under rapid development. Updates and bug fixes are being tracked on the `uga github page`_

.. _`Cython: https://pypi.python.org/pypi
.. _`uga github page`: https://github.com/rmkoesterer/uga

**Current Features**
   - Compatibility with standard `VCFv4.1`_ and `VCFv4.2`_
   - Single SNV association modeling (R base: lm, glm; R `geepack`_: geeglm, R `seqMeta`_: burdenMeta, skatMeta, skatOMeta)
   - Gene/Group based association modeling (with meta analysis)
   - Family based association modeling for single SNV tests
   - Run multiple models as a single submission (alleles are aligned and SNV names need not match)
   - Alignment of compatible SNVs based on position and alleles (A/T and G/C SNVs are ambiguous and are assumed to be pre-aligned)
   - File mapping based on region size or number of SNVs for splitting analyses
   - Automatically split jobs on parallel computing systems using `qsub`_
   - User definable buffered reading for RAM usage control
   - Verification and compilation for parallel distributed jobs
   - `Gzip`_ and `Bgzip`_ / `Tabix`_ mapped output where possible to save disc space
   - A small practice dataset is included with the source code for testing

.. _`VCFv4.1`: http://samtools.github.io/hts-specs/VCFv4.1.pdf
.. _`VCFv4.2`: http://samtools.github.io/hts-specs/VCFv4.2.pdf
.. _`geepack`: https://cran.r-project.org/web/packages/geepack/index.html
.. _`seqMeta`: https://cran.r-project.org/web/packages/seqMeta/index.html
.. _`qsub`: http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html
.. _`Gzip`: http://www.gzip.org/
.. _`Bgzip`: http://www.htslib.org/
.. _`Tabix`: http://www.htslib.org/

**Features Coming Soon**
   - Full documentation
   - Installation with `pip`_
   - `Oxford`_ (3 probabilities for each genotype), `Plink binary`_, and various single allele dosage formatted filetypes
   - Additional association models (R `survival`_: coxph; `lme4`_: lmer, `nlme`_: lme)
   - Family data inclusion in gene/group based tests
   - Interaction terms for relevant models
   - Post modeling meta analysis with genomic control correction
   - Calculation for grouped analysis multiple test correction
   - Publication quality Q-Q and manhattan plots
   - Region-based plots via `Locuszoom`_ software
   - Annotation of results using `SnpEff`_

.. _`Plink binary`: https://www.cog-genomics.org/plink2/input#bed
.. _`Oxford`: http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html
.. _`pip`: https://pypi.python.org/pypi/pip
.. _`survival`: https://cran.r-project.org/web/packages/survival/index.html
.. _`lme4`: https://cran.r-project.org/web/packages/lme4/index.html
.. _`nlme`: https://cran.r-project.org/web/packages/nlme/index.html
.. _`Locuszoom`: http://genome.sph.umich.edu/wiki/LocusZoom_Standalone
.. _`SnpEff`: http://snpeff.sourceforge.net/
.. _`SnpSift`: http://snpeff.sourceforge.net/SnpSift.html

Since parallel computing is sometimes unreliable, analysts need to be able to verify and possibly rerun failed jobs at various stages of the analysis.
In the interest of user efficiency and to avoid limitations induced by excessive automation, uga breaks the analytical process into the following modules.

   - **set** user definable settings
   - **snv** map non-empty regions in genotype/imputed data files
   - **snvgroup** variant and gene/region-based statistical modeling
   - **compile** verify and compile split analysis results
   - **resubmit** resubmit failed jobs for a project
   - **plot** Q-Q and manhattan plots (not yet available)
   - **zoom** region plots (not yet available)
   - **meta** meta-analysis (not yet available)
   - **gc** apply genomic control to 1 or more p-value columns (not yet available)
   - **annot** annotate variant results using `SnpEff`_ and `SnpSift`_ (not yet available)

.. _`SnpEff`: http://snpeff.sourceforge.net/
.. _`SnpSift`: http://snpeff.sourceforge.net/SnpSift.html

Installation
************

This software uses an array of Python modules and R packages. Thus, it may be simpler for users to install it within a clean virtual environment to avoid disrupting system 
Python functionality. The following lists display versions used during development. These modules can be installed easily with `pip`_.

.. _`pip`: https://pypi.python.org/pypi/pip

`Python`_ (2.7.7)

.. _`Python`: https://www.python.org/

Python modules required (may not be part of the Python base install), followed by versions used during development:

   * `singledispatch`_ (3.4.0.3)
   * `rpy2`_ (2.5.2)
   * `multi-key-dict`_ (2.0.1)
   * `numpy`_ (1.9.1)
   * `pandas`_ (0.15.1)
   * `progressbar`_ (2.3)
   * `psutil`_ (2.1.3)
   * `pytabix`_ (0.1)
   * `scipy`_ (0.14.0)
   * `biopython`_ (1.64)
   * `plinkio`_ (0.9.5)
   * `pysam`_ (0.8.2.1)
   * `PyVCF`_ (0.6.7)
   * `Cython`_ (0.22)
   * `XlsxWriter`_ (0.7.2)

.. _`singledispatch`: https://pypi.python.org/pypi/singledispatch
.. _`rpy2`: https://pypi.python.org/pypi/rpy2
.. _`multi-key-dict`: https://pypi.python.org/pypi/multi-key-dict
.. _`numpy`: https://pypi.python.org/pypi/numpy
.. _`pandas`: https://pypi.python.org/pypi/pandas
.. _`progressbar`: https://pypi.python.org/pypi/progressbar
.. _`psutil`: https://pypi.python.org/pypi/psutil
.. _`pytabix`: https://pypi.python.org/pypi/pytabix
.. _`scipy`: https://pypi.python.org/pypi/scipy
.. _`biopython`: https://pypi.python.org/pypi/biopython
.. _`plinkio`: https://pypi.python.org/pypi/plinkio
.. _`pysam`: https://pypi.python.org/pypi/pysam
.. _`PyVCF`: https://pypi.python.org/pypi/PyVCF
.. _`Cython`: https://pypi.python.org/pypi/Cython
.. _`XlsxWriter`: https://pypi.python.org/pypi/XlsxWriter

`R`_ (3.1.1)

.. _`R`: http://www.r-project.org/

R libraries needed for certain analytical and plotting tasks, followed by versions used during development:

   * `ggplot2`_ (1.0.0)
   * `geepack`_ (1.1-6)
   * `lme4`_ (1.1-7)
   * `nlme`_ (1.1-7)
   * `survival`_ (2.37-7)
   * `seqMeta`_ (1.5)
   * `kinship2`_ (1.6.0)

.. _`ggplot2`: http://cran.r-project.org/web/packages/ggplot2/index.html
.. _`geepack`: https://cran.r-project.org/web/packages/geepack/index.html
.. _`seqMeta`: https://cran.r-project.org/web/packages/seqMeta/index.html
.. _`lme4`: http://cran.r-project.org/web/packages/lme4/index.html
.. _`nlme`: https://cran.r-project.org/web/packages/nlme/index.html
.. _`survival`: http://cran.r-project.org/web/packages/survival/index.html
.. _`kinship2`: http://cran.r-project.org/web/packages/kinship2/index.html

Some of these R libraries may have dependencies that need to be installed as well.

Clutter is reduced through consolidation and compression of data and results files via `tabix/bgzip`_ and `gzip`_.

.. _`tabix/bgzip`: http://www.htslib.org/
.. _`gzip`: http://www.gzip.org/

Generating regional plots requires the installation of `locuszoom`_.

.. _`locuszoom`: http://genome.sph.umich.edu/wiki/LocusZoom_Standalone

**Pre-Installation**

To avoid potential errors during installation, you may need to add the location of the R library libR.so file to your BASH_PROFILE 
(ie. .bashrc, .bash_profile, etc). The following command will search your system for this file.
   
   >>> find /usr -name libR.so
	  
Add the resulting path, X, to the following line and add it to your BASH_PROFILE.
   
   export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:X
	  
Make sure you source your BASH_PROFILE again before continuing with the install.
   
   >>> source BASH_PROFILE

**Virtual Environment Preparation**

Installing uga under a Python virtual environment (`virtualenv`_) will ensure that the modules required by uga won't interrupt your system Python install. 
For example, you can install and activate a virtual environment called 'uga-env' as follows:

   >>> mkdir uga-env
   >>> virtualenv -p python uga-env
   >>> source uga-env/bin/activate

.. _`virtualenv`: https://virtualenv.pypa.io/en/latest/

You are now operating a clean base Python installation under a virtual environment.

**Installing uga from source**

Use the following commands to install uga from a source file, uga.tar.gz.

   >>> tar -xvf uga.tar.gz
   >>> cd uga
   >>> pip install -r requirements.txt
   >>> python setup.py install

**Installing uga with pip (not yet available)**

The simplest way to install uga is with `pip`_, as follows.

   >>> pip install uga

.. _`pip`: https://pypi.python.org/pypi/pip

**Note**: If you install uga under a virtual environment, you need to source the environment as shown above before running any task in uga.

   >>> source uga-env/bin/activate

Verify that uga is functional using the following command to display help.

   >>> uga -h

**Parallel computing**

While you may simply run uga on a single cpu system, if you have access to a parallel computing cluster or even a single multiple core
processor, you will be able to take advantage of the self-managed parallel mode of use for which this software was designed. 
This release was tested on a system which deploys Sun Grid Engine and `qsub`_ for job management and will likely be compatible 
with other PBS systems.

.. _`qsub`: http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

References
==========

Manuscript to be submitted

Contact
=======

- **Author**: `Ryan Koesterer`_

.. _`Ryan Koesterer`: https://github.com/rmkoesterer/uga

License
=======

Universal Genome Analyst (uga) is distributed under the GNU General Public License v3:
   
   Copyright (c) 2015 Ryan Koesterer

   This program is free software: you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation, either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see
   <http://www.gnu.org/licenses/>
