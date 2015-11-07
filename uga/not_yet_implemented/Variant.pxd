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
cimport numpy as np
cimport cython
import math

cdef double CalcCallrate(np.ndarray[np.float64_t, ndim=1])
cdef double CalcFreq(np.ndarray[np.float64_t, ndim=1])
cdef double CalcMAC(np.ndarray[np.float64_t, ndim=1])
cdef double CalcFreqX(np.ndarray[np.float64_t, ndim=1] male, np.ndarray[np.float64_t, ndim=1] female)
cdef double CalcRsq(np.ndarray[np.float64_t, ndim=1])
