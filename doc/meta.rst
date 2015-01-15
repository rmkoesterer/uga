Meta-Analysis
*************

	The meta module performs meta-analysis and requires a configuration file.

.. argparse::
   :ref: uga.Parse.Parser
   :prog: uga
   :path: meta

Configuration Files
===================

	The meta module requires a configuration file to define parameters for the meta-analysis. 
	
Required
--------
	
	out phenotype.meta
		
Column Names
------------
	
	chr chr

	pos pos

	marker marker

	a1 a1

	a2 a2

	freq freq

	effect marker.effect

	stderr marker.stderr

	or marker.or

	z marker.z

	p marker.p

Filters
-------

	filter freq >= 0.03
	filter freq <= 0.97

Tags
----

	tag aa

Sample Size
-----------

	n 2180

Genomic Control
---------------
	
	gc 1.232

Process File Name
-----------------

	process_file assoc.aa.chr1bp900000-1000000.gz
