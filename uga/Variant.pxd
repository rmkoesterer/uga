## Copyright (c) 2015 Ryan Koesterer GNU General Public License v3
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import pandas as pd
import Geno
cimport Geno
cimport numpy as np
cimport cython
import math

cdef class Ref:
	cdef public object db
	cpdef load(self, np.ndarray v)
	cpdef update(self, np.ndarray v)
cpdef complement(allele)
cpdef get_universal_variant_id(chr_py,pos_py,a1_py,a2_py,delim_py)
cdef double calc_callrate(np.ndarray[np.float64_t, ndim=1])
cdef double calc_freq(np.ndarray[np.float64_t, ndim=1])
cdef double calc_mac(np.ndarray[np.float64_t, ndim=1])
cdef double calc_freqx(np.ndarray[np.float64_t, ndim=1] male, np.ndarray[np.float64_t, ndim=1] female)
cdef double calc_rsq(np.ndarray[np.float64_t, ndim=1])
