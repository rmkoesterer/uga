Universal Genome Analyst
========================
  
Universal Genome Analyst (uga) is a toolbox designed to assist biomedical researchers in complex genomic data analysis, making use of many powerful existing 
R packages and Python modules along with large scale computing cluster integration to provide the following features.

* Compatibility with standard VCF4.0/4.1, Oxford (Impute2 output with 3 probabilities for each genotype), Plink binary, and various single allele dosage format files
* Mapping of files based on region size or number of markers for splitting analyses
* Automatic deployment of jobs on parallel computing systems using qsub
* Verification and compilation of parallel distributed jobs
* Marker association modelling:
   generalized linear models
   linear mixed effects models
   generalized estimating equations
   survival analysis
* Gene/Locus based association modelling
   effective test correction
   burden test
   sequence kernel association test (SKAT)
   optimal unified sequence kernel association test (SKAT-O)
   family Based SKAT and SKAT-O
* Meta analysis
   uses unique marker naming based on both position and alleles to allow compatibility between multiple marker naming conventions
* automatically aligns compatible markers alleles, eliminating the need to align input results manually
* Optional genomic inflation correction (genomic control)
* Publication quality Q-Q and Manhattan Plots
* Results filtering and reporting
* Gzip and Bgzip / Tabix mapped output where possible to save disc space
* User definable buffered reading for RAM usage control
* Report compatibility for annotation with SnpEff and SnpSift

Installation
============

This software connects a vast array of Python and R based packages that may lead to incompatibilities in your system Python installation. Thus, it might be simpler for users
to install a Python virtual environment to avoid disrupting system Python functionality. Additionally, uga requires an installation of R and a few packages used for analysis 
and plotting. The following lists are necessary tools that uga needs to perform all tasks. Python modules are easiest installed using pip as described in the section labeled 
'Installation', thus you may skip ahead to the Installation section to begin. The following lists display versions used during development.

Python (2.7.7)

Python modules required (may not be part of the Python base install), followed by versions used during development:

* singledispatch (3.4.0.3)
* rpy2 (2.5.2)
* multi-key-dict (2.0.1)
* numpy (1.9.1)
* pandas (0.15.1)
* progressbar (2.3)
* psutil (2.1.3)
* pytabix (0.1)
* scipy (0.14.0)
* biopython (1.64)
* plinkio (0.9.5)

R (3.1.1)

R libraries needed for certain analytical and plotting tasks, followed by versions used during development:

* ggplot2 (1.0.0)
* geepack (1.1-6)
* lme4 (1.1-7)
* survival (2.37-7)
* seqMeta (1.5)
  
Some of these libraries may have dependencies which will need to be installed as well.
   
Tabix and bgzip

In order to reduce file clutter and encourage the consolidation and compression of data and results files, uga makes extensive use of both Tabix and bgzip. 
These tools are generally release as part of the `Samtools`_ suite
	
.. _`Samtools`: http://www.htslib.org/

**Installation Steps**

Before installing uga, for full functionality, the following should be installed and working (see above for development version information):

* Python (base install)

* virtualenv

* R (base install), plus the following libraries
   ggplot2
   geepack
   lme4
   survival
   seqMeta
* Tabix
* bgzip
   
In order to avoid potential errors during installation, you may need to add the location of the R library libR.so file to your BASH_PROFILE 
(ie. .bashrc, .bash_profile, etc). The following command will search your system for this file.
   
   >>> find /usr -name libR.so
	  
Add the resulting path, X, the following line and add the line to your BASH_PROFILE.
   
   export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:X
	  
Make sure you source your BASH_PROFILE again before continuing with the install.
   
   >>> source BASH_PROFILE
	  
Choose a directory in which you'd like your virtual environment to live, for example 'uga-env', then install the environment and source it.

   >>> mkdir uga-env
   >>> virtualenv -p python uga-env
   >>> source uga-env/bin/activate
  
After sourcing your virtual environment, you can install the required Python modules for uga as follows.

   >>> cd uga
   >>> pip install -r requirements.txt
   >>> pip install uga

* There is a qsub wrapper included in your installation directory (bin/.uga_wrapper.py). This needs to be copied or moved to your home directory to allow uga to submit
jobs to your computing cluster using the qsub command.
	  
Note: The virtual environment created during installation is the environment under which uga must be run, thus you need to source the environment
before running any task in uga.

**Parallel computing**

While you may simply run uga on a single cpu system, if you have access to a parallel computing cluster, 
you will be able to take advantage of the self-managed parallel mode of use for which this software was designed. 
This release was tested on a system which deploys Sun Grid Engine for job management, but simple modifications to
the uga_submit.py script may allow the use of other PBS systems, such as Torque.

References
==========

Manuscript to be submitted

Contact
=======

- **Author**: `Ryan Koesterer`_

`Documentation`_

.. _`Ryan Koesterer`: uga-feedback@gmail.com
.. _`Documentation`: http://rmkoesterer.github.io/uga-doc/
