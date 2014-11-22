Universal Genome Analyst
========================
  
- **Authors**: `Ryan Koesterer`_

.. _`Ryan Koesterer`: koesterr@bu.edu
.. _`Boston University Biomedical Genetics`: http://www.bumc.bu.edu/genetics

Features
--------

- Marker Analysis (implementation in R)
   - generalized linear models (glm)
   - linear mixed effects models (lme)   
   - generalized estimating equations (geeglm)
   - Cox proportional hazards (coxph)
- Meta analysis
- Annotation
- Plotting


Dependencies
------------
 To check for these required python modules, type 'pydoc modules' on the commandline and search for them by name in the resulting list 
 - SGE parallel computing environment
- R 2.15.3
- python2.7
easy_install virtualenv --user
cd into * (some directory where you'd like your virtual environment to live)
virtualenv * (if currently using 2.7)
	virtualenv -p python2.7 * (if not currently using 2.7)
source */bin/activate
	pip install pytabix
	pip install rpy2
	pip install multi_key_dict
	pip install numpy
	pip install progressbar
	pip install pandas
	pip install psutil
	
- reshape (R package)
- grDevices (R package)
- gtools (R package)
- geepack (R package)
- lme4 (R package)
- survival (R package)

Use the following command to display library links for your virtual python environment
ldd */bin/python
make sure the value, X, in the line 'libpython2.7.so.1.0 ==> X' is added to LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:X

Python module installation
--------------------------

Use the following method to install python modules::

 python2.7 PATH_TO_MODULE/setup.py install --user
 
 find /usr -name libR.so
 
 export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:PATH_TO_DIR_CONTAINING_libR.so

Additions to ~/.bashrc
----------------------

 module load R/R-2.15.3_gnu-4.4.6

 module load tabix/0.2.6
 
 module load python2.7/Python-2.7.3_gnu446

 export TEMPDIR=/scratch
 
 export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/IT/R-2.15.3/lib64/R/lib/
