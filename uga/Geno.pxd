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
import Variant
cimport Variant
cimport numpy as np
cimport cython
import math

cdef class Variants:
	cdef public bytes filename, sample_filename, region, group_id, method
	cdef public object handle, region_iter, snvgroup_map, results, header, cols, dtypes, snv_results_tagged, duplicated
	cdef public unsigned int chr, start, end
	cdef public np.ndarray genos, data, info, snv_chunk, snvgroup_chunk, snv_results, samples
	cpdef align(self, Variant.Ref ref)
