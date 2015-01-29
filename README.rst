Universal Genome Analyst
========================
  
- **Authors**: `Ryan Koesterer`_

.. _`Ryan Koesterer`: koesterr@bu.edu
.. _`Boston University Biomedical Genetics`: http://www.bumc.bu.edu/genetics

Features
--------

- Marker Analysis (R packages)
   - generalized linear models (glm)
   - linear mixed effects models (lme)   
   - generalized estimating equations (geeglm)
   - Cox proportional hazards (coxph)
- Gene Based Analysis
   - Li and Ji gene based test (top snp in each gene corrected for number of effective tests)
- Meta analysis


Dependencies
------------

 - SGE or similar parallel computing environment
 
 - R 3.1.1

 - python2.7
	
 - reshape (R package)

 - grDevices (R package)

 - gtools (R package)

 - geepack (R package)

 - lme4 (R package)

 - survival (R package)

Installation
------------
 
 easy_install virtualenv --user

 cd into \* (some directory where you'd like your virtual environment to live)

 virtualenv -p python2.7 \*

 source \*/bin/activate

 pip install -r requirements.txt

 python setup.py install

# 
cd uga-2015.01a
mkdir env
module load python_modules/virtualenv
virtualenv -p python env
source env/bin/activate
cp ../uga-2015.01/* .
cp -r ../uga-2015.01/uga .
cp -r ../uga-2015.01/bin .
pip install -r requirements.txt
python setup.py install
#

Hints
-----

If you receive an error regarding libR.so, use the following command to find it's location

 find /usr -name libR.so
 
Use the following command to display library links for your virtual python environment

 ldd \*/bin/python

make sure the value, X, in the line 'libpython2.7.so.1.0 ==> X' is added to LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:X in your bashrc
 
Suggested additions to ~/.bashrc
--------------------------------

Add these to your bashrc file
 
 export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:X
 
where X is the path to the directory containing libR.so (see above)

 export TEMPDIR=/scratch
 
 export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:PATH_TO_R_LIB
