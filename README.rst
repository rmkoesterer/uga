Universal Genome Analyst
************************

Universal Genome Analyst (**uga**) is an open, flexible, and efficient tool for the distribution, management, and visualization of whole genome data analyses. 
It is designed to assist biomedical researchers in complex genomic data analysis through the use of a low level interface between the powerful R statistical 
environment and Python, allowing for rapid integration of emerging analytical strategies. This project uses `Cython`_ for a significant reduction in computation 
time and researchers with access to a high performance computing cluster or with access to multiple cores will find time-saving features for parallel analysis 
using a flexible, yet controlled, commandline interface.

.. _`Cython: https://pypi.python.org/pypi

This software is currently under rapid development. Updates and bug fixes are being tracked on the `uga github page`_

.. _`uga github page`: https://github.com/rmkoesterer/uga

**Current Features**
   - Compatibility with standard `VCFv4.1`_ and `VCFv4.2`_ files
   - Single variant association modeling (R base: lm, glm; R `geepack`_: geeglm, R `seqMeta`_: singlesnpMeta)
   - Gene/Group based association modeling (with meta analysis: R `seqMeta`_: burdenMeta, skatMeta, skatOMeta)
   - Family based association modeling
   - Run multiple models as a single submission (variant names need not match)
   - Publication quality Q-Q and manhattan plots
   - Alignment of compatible variants based on genomic position and both alleles (A/T and G/C SNVs are ambiguous and are assumed to be pre-aligned)
   - An array of job splitting capabilities
   - Automatically submit jobs on parallel computing systems using `qsub`_
   - multiple processor parallelization in addition to cluster parallelization
   - User definable buffered reading for RAM usage control
   - Verification and compilation for parallel distributed jobs
   - `Gzip`_ and `Bgzip`_ / `Tabix`_ mapped output where possible to save disc space

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
   - Additional association models (R `survival`_: coxph; `lme4`_: lmer, `nlme`_: lme)
   - Family data inclusion in gene/group based tests
   - Genomic control correction
   - Post modeling meta analysis
   - Calculation for grouped analysis multiple test correction
   - Region-based plots via `Locuszoom`_ software
   - Results annotation using `SnpEff`_

.. _`survival`: https://cran.r-project.org/web/packages/survival/index.html
.. _`lme4`: https://cran.r-project.org/web/packages/lme4/index.html
.. _`nlme`: https://cran.r-project.org/web/packages/nlme/index.html
.. _`Locuszoom`: http://genome.sph.umich.edu/wiki/LocusZoom_Standalone
.. _`SnpEff`: http://snpeff.sourceforge.net/

Since parallel computing is sometimes unreliable, analysts need to be able to verify and possibly rerun failed jobs at various stages of the analysis.
In the interest of user efficiency and to avoid limitations induced by excessive automation, uga breaks the analytical process into the following modules.

   - **set** user definable settings
   - **snv** single variant statistical modeling
   - **snvgroup** gene/region-based statistical modeling
   - **compile** verify and compile split analysis results
   - **resubmit** automatically resubmit failed jobs for a project
   - **snvplot** Q-Q and manhattan plots for single snv tests
   - **snvgroupplot** Q-Q and manhattan plots for snv group tests
   - **zoom** regional plots (not yet available)
   - **meta** meta-analysis (not yet available)
   - **gc** apply genomic control to results (not yet available)
   - **annot** annotate variant results using `SnpEff`_ (not yet available)

.. _`SnpEff`: http://snpeff.sourceforge.net/

Installation
************

This software uses a variety of Python modules, R packages, and some stand-alone software. Thus, the easiest method for installation is using `conda`_.
The required modules are installed easily using the environment.yml file included with this distribution, as described in the Pre-installation section below.

.. _`conda`: http://conda.pydata.org/docs/

Consolidation and compression of data and results files requires `tabix/bgzip`_ and `gzip`_.

.. _`tabix/bgzip`: http://www.htslib.org/
.. _`gzip`: http://www.gzip.org/

Generating regional plots requires the installation of `locuszoom`_.

.. _`locuszoom`: http://genome.sph.umich.edu/wiki/LocusZoom_Standalone

Annotation of results requires `SnpEff`_.

.. _`SnpEff`: http://snpeff.sourceforge.net/

**Pre-Installation (preparing the environment for uga)**

To prepare your system for uga, you need to `clone an environment`_ using conda. You will need the included environment.yml file from the source code and a number of 
conda packages from `my anaconda cloud channel`_ and another custom channel (both are listed in the environment.yml file). After downloading the most recent 
release (available `here`_), use the following commands to begin the installation.

.. _`clone an environment`: http://conda.pydata.org/docs/using/envs.html#clone-an-environment
.. _`my anaconda cloud channel`: https://conda.anaconda.org/rmkoesterer
.. _`here`: https://github.com/rmkoesterer/uga/releases

   >>> tar -xvf uga.tar.gz
   >>> cd uga
   >>> conda env create -f environment.yml
   >>> source activate uga_python2.7

**Installing uga from source**

Once you have successfully prepared the conda environment, use the following command to install uga from source.

   >>> python setup.py install

**Cutting Edge Install**

Keeping up with the most current changes may be of interest to you as I will be rapidly adding features well into 2016. Thus, under the realization 
that you may encounter unexpected behavior and bugs, you may want to run a fork of this repository rather than installing from source. See this tutorial describing
how to `fork this repository`_

.. _`fork this repository`: https://help.github.com/articles/fork-a-repo/

**Note**: If you install uga under a conda environment, you need to source the environment as shown above before running any task in uga.

   >>> source activate uga_python2.7

Verify that uga is functional using the following command to display help.

   >>> uga -h

Note: further help is provided after selecting a specific module, ie.

   >>> uga snv -h

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

Please report any bugs or issues using the `Issues`_ tab on this page. I will respond to all concerns as quickly as possible.

.. _`Issues`: https://github.com/rmkoesterer/uga/issues

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
