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

!> Definition of the abstract base dispersion model for the evaluation of C6 coefficients.
module dftd4_model_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use multicharge, only : mchrg_model_type
   use dftd4_cache, only : dispersion_cache
   implicit none
   private

   public :: dispersion_model, d4_qmod


   !> Abstract base dispersion model to evaluate C6 coefficients
   type, abstract :: dispersion_model

      !> Number of atoms coupled to by pairwise parameters
      integer :: ncoup

      !> Number of frequency grid points for dynamic polarizabilities
      integer :: ngrid

      !> Charge scaling height
      real(wp) :: ga

      !> Charge scaling steepness
      real(wp) :: gc

      !> Effective nuclear charges
      real(wp), allocatable :: zeff(:)

      !> Chemical hardness
      real(wp), allocatable :: eta(:)

      !> Electronegativity
      real(wp), allocatable :: en(:)

      !> Covalent radii for coordination number
      real(wp), allocatable :: rcov(:)

      !> Expectation values for C8 extrapolation
      real(wp), allocatable :: r4r2(:)

      !> Number of reference systems
      integer, allocatable :: ref(:)

      !> Number of Gaussian weights for each reference
      integer, allocatable :: ngw(:, :)

      !> Reference coordination numbers
      real(wp), allocatable :: cn(:, :)

      !> Reference partial charges
      real(wp), allocatable :: q(:, :)

      !> Reference dynamic polarizabilities
      real(wp), allocatable :: aiw(:, :, :)

      !> Reference C6 coefficients
      real(wp), allocatable :: c6(:, :, :, :)

      !> Multicharge model
      class(mchrg_model_type), allocatable :: mchrg 

   contains

      !> Update cache with dispersion coefficients and properties
      procedure(update), deferred :: update

      !> Evaluate atomic polarizabilities from cache
      procedure(get_polarizabilities), deferred :: get_polarizabilities

   end type dispersion_model

   abstract interface

      !> Update dispersion cache with precomputed coefficients and properties
      subroutine update(self, mol, cache, cn, q, grad, only_c6)
         import dispersion_model, dispersion_cache, structure_type, wp
         !> Instance of the dispersion model
         class(dispersion_model), intent(in) :: self
         !> Molecular structure data
         class(structure_type), intent(in) :: mol
         !> Dispersion cache to populate
         type(dispersion_cache), intent(inout) :: cache
         !> Coordination number of every atom
         real(wp), intent(in) :: cn(:)
         !> Partial charge of every atom
         real(wp), intent(in) :: q(:)
         !> Whether to compute derivatives
         logical, intent(in), optional :: grad
         !> Whether to compute only C6 coefficients
         logical, intent(in), optional :: only_c6
      end subroutine update

      !> Calculate atomic polarizabilities from cache
      subroutine get_polarizabilities(self, cache, alpha, alphaqq, &
         & dadcn, dadq, daqqdcn, daqqdq)
         import dispersion_model, dispersion_cache, wp
         !> Instance of the dispersion model
         class(dispersion_model), intent(in) :: self
         !> Dispersion cache containing polarizabilities
         type(dispersion_cache), intent(in) :: cache
         !> Static dipole-dipole polarizabilities for all atoms
         real(wp), intent(out) :: alpha(:)
         !> Static quadrupole-quadrupole polarizabilities for all atoms
         real(wp), intent(out) :: alphaqq(:)
         !> Derivative of dipole polarizibility w.r.t. coordination number
         real(wp), intent(out), optional :: dadcn(:)
         !> Derivative of dipole polarizibility w.r.t. partial charge
         real(wp), intent(out), optional :: dadq(:)
         !> Derivative of quadrupole polarizibility w.r.t. coordination number
         real(wp), intent(out), optional :: daqqdcn(:)
         !> Derivative of quadrupole polarizibility w.r.t. partial charge
         real(wp), intent(out), optional :: daqqdq(:)
      end subroutine get_polarizabilities

   end interface


   !> Possible reference charges for D4
   type :: enum_qmod

      !> Electronegativity equilibration charges
      integer :: eeq = 1

      !> GFN2-xTB Mulliken partial charges
      integer :: gfn2 = 2

      !> Bond-Capcity Electronegativity equilibration charges
      integer :: eeqbc = 3

   end type enum_qmod

   !> Actual enumerator for D4 reference charges
   type(enum_qmod), parameter :: d4_qmod = enum_qmod()
   !DEC$ ATTRIBUTES DLLEXPORT :: d4_qmod

end module dftd4_model_type
