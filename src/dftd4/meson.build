# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

subdir('damping')
subdir('data')
subdir('model')

srcs += files(
  'blas.F90',
  'charge.f90',
  'cutoff.f90',
  'damping.f90',
  'data.f90',
  'disp.f90',
  'model.f90',
  'ncoord.f90',
  'numdiff.f90',
  'output.f90',
  'param.f90',
  'reference.f90',
  'utils.f90',
  'version.f90',
)
if get_option('api') or get_option('python')
  srcs += files('api.f90')
endif
if get_option('api_v2')
  srcs += files('compat.f90')
endif
