! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Generic interface to define damping functions for the DFT-D models
module dftd4_damping
   use mctc_env, only : wp
   implicit none
   private

   public :: damping_param

   !> Abstract base type for damping parameterizations (undamped case)
   type, abstract :: damping_param
      !> Scaling factor for C6/R^6 term
      real(wp) :: s6 = 1.0_wp
      !> Scaling factor for C8/R^8 term
      real(wp) :: s8
      !> Scaling factor for C9/R^9 term
      real(wp) :: s9 = 1.0_wp
      !> Linear damping radius dependence
      real(wp) :: a1
      !> Constant damping radius shift
      real(wp) :: a2
      !> Zero-damping parameter
      real(wp) :: alp = 16.0_wp
   contains
      !> Evaluate two-body damping factor
      procedure :: get_2b_damp
      !> Evaluate two-body damping factor with derivarives
      procedure :: get_2b_derivs
      !> Evaluate three-body damping factor
      procedure :: get_3b_damp
      !> Evaluate three-body damping factor with derivarives
      procedure :: get_3b_derivs
   end type damping_param

contains

!> Evaluate two-body damping factor
pure subroutine get_2b_damp(self, r2, rdamp, d6, d8)
   !> Damping parameters
   class(damping_param), intent(in) :: self
   !> Square of interatomic distance
   real(wp), intent(in) :: r2
   !> Damping radius
   real(wp), intent(in) :: rdamp
   !> Damping factor for C6/R^6 term
   real(wp), intent(out) :: d6
   !> Damping factor for C8/R^8 term
   real(wp), intent(out) :: d8

   d6 = 1.0_wp / r2**3
   d8 = 1.0_wp / r2**4

end subroutine get_2b_damp

!> Evaluate two-body damping factor with derivarives
pure subroutine get_2b_derivs(self, r2, rdamp, d6, d8, d6dr, d8dr)
   !> Damping parameters
   class(damping_param), intent(in) :: self
   !> Square of interatomic distance
   real(wp), intent(in) :: r2
   !> Damping radius
   real(wp), intent(in) :: rdamp
   !> Damping factor for C6/R^6 term
   real(wp), intent(out) :: d6
   !> Damping factor for C8/R^8 term
   real(wp), intent(out) :: d8
   !> Derivative of damping factor for C6/R^6 w.r.t the interatomic distance
   real(wp), intent(out) :: d6dr
   !> Derivative of damping factor for C8/R^8 w.r.t the interatomic distance
   real(wp), intent(out) :: d8dr

   d6 = 1.0_wp / r2**3
   d8 = 1.0_wp / r2**4

   d6dr = -6*r2**2*d6**2
   d8dr = -8*r2**3*d8**2

end subroutine get_2b_derivs

!> Evaluate three-body damping factor
pure subroutine get_3b_damp(self, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9)
   !> Rational damping parameters
   class(damping_param), intent(in) :: self
   !> Product of pairwise interatomic distances
   real(wp), intent(in) :: r
   !> Square of interatomic distance between atoms i and j
   real(wp), intent(in) :: r2ij
   !> Square of interatomic distance between atoms i and k
   real(wp), intent(in) :: r2ik
   !> Square of interatomic distance between atoms j and k
   real(wp), intent(in) :: r2jk
   !> Product of pairwise damping radii
   real(wp), intent(in) :: rdamp
   !> Pairwise damping radius of atoms i and j
   real(wp), intent(in) :: rdampij
   !> Pairwise damping radius of atoms i and k
   real(wp), intent(in) :: rdampik
   !> Pairwise damping radius of atoms j and k
   real(wp), intent(in) :: rdampjk
   !> Damping factor for C9/R^9 term
   real(wp), intent(out) :: d9

   d9 = 1.0_wp / r**3

end subroutine get_3b_damp

!> Evaluate three-body damping factor
pure subroutine get_3b_derivs(self, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9, d9drij, d9drik, d9drjk)
   !> Rational damping parameters
   class(damping_param), intent(in) :: self
   !> Product of pairwise interatomic distances
   real(wp), intent(in) :: r
   !> Square of interatomic distance between atoms i and j
   real(wp), intent(in) :: r2ij
   !> Square of interatomic distance between atoms i and k
   real(wp), intent(in) :: r2ik
   !> Square of interatomic distance between atoms j and k
   real(wp), intent(in) :: r2jk
   !> Product of pairwise damping radii
   real(wp), intent(in) :: rdamp
   !> Pairwise damping radius of atoms i and j
   real(wp), intent(in) :: rdampij
   !> Pairwise damping radius of atoms i and k
   real(wp), intent(in) :: rdampik
   !> Pairwise damping radius of atoms j and k
   real(wp), intent(in) :: rdampjk
   !> Damping factor for C9/R^9 term
   real(wp), intent(out) :: d9
   !> Derivative of damping factor for C9/R^9 
   !> w.r.t the interatomic distance between atoms i and j
   real(wp), intent(out) :: d9drij
   !> Derivative of damping factor for C9/R^9 
   !> w.r.t the interatomic distance between atoms i and k
   real(wp), intent(out) :: d9drik
   !> Derivative of damping factor for C9/R^9 
   !> w.r.t the interatomic distance between atoms j and k
   real(wp), intent(out) :: d9drjk

   real(wp) :: dftmp

   d9 = 1.0_wp / r**3

   d9drij = -3 * d9 / r2ij
   d9drik = -3 * d9 / r2ik
   d9drjk = -3 * d9 / r2jk

end subroutine get_3b_derivs

end module dftd4_damping
