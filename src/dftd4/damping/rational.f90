! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
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

!> Implementation of the rational (Becke--Johnson) damping function.
module dftd4_damping_rational
   use dftd4_damping, only : damping_param
   use mctc_env, only : wp
   implicit none
   private

   public :: rational_damping_param


   !> Rational (Becke-Johnson) damping model
   type, extends(damping_param) :: rational_damping_param
   contains
      !> Evaluate two-body damping factor
      procedure :: get_2b_damp
      !> Evaluate two-body damping factor with derivarives
      procedure :: get_2b_derivs
      !> Evaluate three-body damping factor
      procedure :: get_3b_damp
      !> Evaluate three-body damping factor with derivarives
      procedure :: get_3b_derivs
   end type rational_damping_param


contains


!> Evaluate two-body damping factor
pure subroutine get_2b_damp(self, r2, rdamp, d6, d8)
   !> Rational damping parameters
   class(rational_damping_param), intent(in) :: self
   !> Square of interatomic distance
   real(wp), intent(in) :: r2
   !> Damping radius
   real(wp), intent(in) :: rdamp
   !> Damping factor for C6/R^6 term
   real(wp), intent(out) :: d6
   !> Damping factor for C8/R^8 term
   real(wp), intent(out) :: d8

   d6 = self%s6 / (r2**3 + rdamp**6)
   d8 = self%s8 / (r2**4 + rdamp**8)

end subroutine get_2b_damp


!> Evaluate two-body damping factor with derivarives
pure subroutine get_2b_derivs(self, r2, rdamp, d6, d8, d6dr, d8dr)
   !> Rational damping parameters
   class(rational_damping_param), intent(in) :: self
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

   d6 = self%s6 / (r2**3 + rdamp**6)
   d8 = self%s8 / (r2**4 + rdamp**8)

   d6dr = -6.0_wp * r2**2 * d6**2 / self%s6
   d8dr = -8.0_wp * r2**3 * d8**2 / self%s8

end subroutine get_2b_derivs


!> Evaluate three-body damping factor
pure subroutine get_3b_damp(self, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9)
   !> Rational damping parameters
   class(rational_damping_param), intent(in) :: self
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

   d9 = self%s9 / (1.0_wp + 6.0_wp * (rdamp / r)**(self%alp / 3.0_wp))

end subroutine get_3b_damp


!> Evaluate three-body damping factor
pure subroutine get_3b_derivs(self, r, r2ij, r2ik, r2jk, &
   & rdamp, rdampij, rdampik, rdampjk, d9, d9drij, d9drik, d9drjk)
   !> Rational damping parameters
   class(rational_damping_param), intent(in) :: self
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

   d9 = self%s9 / (1.0_wp + 6.0_wp * (rdamp / r)**(self%alp / 3.0_wp))

   dftmp = -2.0_wp * self%alp * (rdamp / r)**(self%alp / 3.0_wp) * d9**2
   d9drij = dftmp / r2ij
   d9drik = dftmp / r2ik
   d9drjk = dftmp / r2jk

end subroutine get_3b_derivs

end module dftd4_damping_rational
