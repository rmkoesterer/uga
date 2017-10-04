Universal Genome Analyst
************************

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.578712.svg
   :target: https://doi.org/10.5281/zenodo.578712

Universal Genome Analyst (**uga**) is an open, flexible, and efficient tool for the distribution, management, and visualization of whole genome data analyses. 
It is designed to assist biomedical researchers in complex genomic data analysis through the use of a low level interface between the powerful R statistical 
environment and Python, allowing for rapid integration of emerging analytical strategies. This project uses `Cython`_ for a significant reduction in computation 
time and researchers with access to a high performance computing cluster or with access to multiple cores will find time-saving features for parallel analysis 
using a flexible, yet controlled, commandline interface.

.. _`Cython`: https://pypi.python.org/pypi

This software is currently under rapid development. Updates and bug fixes are being tracked on the `uga github page`_

.. _`uga github page`: https://github.com/rmkoesterer/uga

**Notable Features**
   - Single variant association modeling (R base: lm, glm; R `geepack`_: geeglm, R `seqMeta`_: singlesnpMeta, R `lme4`_: lmer)
   - Gene/Group based association modeling (with meta analysis: R `seqMeta`_: burdenMeta, skatMeta, skatOMeta)
   - Family based single variant association modeling
   - Publication quality Q-Q and manhattan plots
   - Genomic control correction
   - Post modeling meta analysis
   - Run multiple models as a single submission (variant names need not match)
   - Alignment of compatible variants based on genomic position and both alleles (A/T and G/C SNVs are ambiguous and are assumed to be pre-aligned)
   - Automatic job splitting (with job array queueing)
   - Input data split by chromosome can be linked via wildcard
   - Automatically submit jobs on parallel computing systems using `qsub`_
   - multiple processor parallelization in addition to cluster parallelization
   - User definable buffered reading for RAM usage control
   - Verification and compilation for parallel distributed jobs
   - `Gzip`_ and `Bgzip`_ / `Tabix`_ mapped output where possible to save disc space

.. _`geepack`: https://cran.r-project.org/web/packages/geepack/index.html
.. _`seqMeta`: https://cran.r-project.org/web/packages/seqMeta/index.html
.. _`lme4`: https://cran.r-project.org/web/packages/lme4/index.html
.. _`qsub`: http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html
.. _`Gzip`: http://www.gzip.org/
.. _`Bgzip`: http://www.htslib.org/
.. _`Tabix`: http://www.htslib.org/

**Planned For Future Releases**
   - Full documentation
   - Additional association models (R `survival`_: coxph; `nlme`_: lme)
   - Family data inclusion in gene/group based tests
   - Calculation for grouped analysis multiple test correction
   - Region-based plots via `Locuszoom`_ software
   - Results annotation using `SnpEff`_

.. _`survival`: https://cran.r-project.org/web/packages/survival/index.html
.. _`nlme`: https://cran.r-project.org/web/packages/nlme/index.html
.. _`Locuszoom`: http://genome.sph.umich.edu/wiki/LocusZoom_Standalone
.. _`SnpEff`: http://snpeff.sourceforge.net/

Since parallel computing is sometimes unreliable, analysts need to be able to verify and possibly rerun failed jobs at various stages of the analysis.
In the interest of user efficiency and to avoid limitations induced by excessive automation, uga breaks the analytical process into the following modules.

   - **settings** user definable settings
   - **snv** single variant statistical modeling
   - **snvgroup** gene/region-based statistical modeling
   - **meta** meta-analysis
   - **compile** verify and compile split analysis results
   - **resubmit** automatically resubmit failed jobs for a project
   - **snvplot** Q-Q and manhattan plots for snv tests
   - **snvgroupplot** Q-Q and manhattan plots for snvgroup tests
   - **filter** filter results / apply genomic control to results
   - **merge** merge and annotate results with external files
   - **tools** run any command line tool with ability to include genomic region automatically

Installation
************

This software uses a variety of Python modules, R packages, and some stand-alone software. Thus, the easiest method for installation is to use one of two platforms of the 
software `conda`_; either `Anaconda`_ or `Miniconda`_.

.. _`conda`: https://conda.io/docs/download.html
.. _`Anaconda`: https://www.continuum.io/downloads
.. _`Miniconda`: https://conda.io/miniconda.html

Also, consolidation and compression of data and results files requires `tabix/bgzip`_ and `gzip`_.

.. _`tabix/bgzip`: http://www.htslib.org/
.. _`gzip`: http://www.gzip.org/

To prepare your system for uga, you need to `clone an environment`_. You will need the included environment.yml file from the source code and a number of 
packages from `my anaconda cloud channel`_ and other custom channels (listed in the environment.yml file). After downloading the most recent 
release (available `here`_), use the following commands to begin the installation.

.. _`clone an environment`: http://conda.pydata.org/docs/using/envs.html#clone-an-environment
.. _`my anaconda cloud channel`: https://conda.anaconda.org/rmkoesterer
.. _`here`: https://github.com/rmkoesterer/uga/releases

For the sake of this tutorial, let's assume the release version is 'X'.

   >>> tar -xvf uga-X.tar.gz
   >>> cd uga-X

At this point you may change the name of the environment to anything you'd prefer by modifying the first line of the environment.yml file. For these instructions, we will 
assume the name is unchanged from 'uga'.

   >>> conda env create -f environment.yml
   >>> source activate uga

Now that your environment is activated, you are ready to install uga from source.

   >>> python setup.py install

**Cutting Edge Install**

Keeping up with the most current changes may be of interest to you as I will likely continue to add features and fix bugs on a regular basis. Thus, you may want to run a fork 
of this repository rather than installing from source. See a tutorial describing how to `fork this repository`_.

.. _`fork this repository`: https://help.github.com/articles/fork-a-repo/

Getting Started
***************

If you install uga under a conda environment, you need to source the environment as shown above before running any task in uga.

   >>> source activate uga

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

Citation
========

Please cite this software as follows. A manuscript is in the works and yet to be submitted.

Koesterer, Ryan. Universal Genome Analyst (uga). https://github.com/rmkoesterer/uga. DOI: 10.5281/zenodo.578712.

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
