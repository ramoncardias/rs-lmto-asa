!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!   
! MODULE: Self consistent
!   
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!   
! DESCRIPTION: 
!> Module to handle data relative to basic control operations
!------------------------------------------------------------------------------

module self_consistent_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use precision_mod
  use string_mod
  implicit none

  private

  !> Module's main structure
  type, public :: self_consistent

  contains
    procedure :: run
    final :: destructor
  end type self_consistent


contains

  subroutine destructor(this)
    type(self_consistent) :: this
  end subroutine destructor

  subroutine run(this)
    class(self_consistent) :: this
  end subroutine run
    
end module self_consistent_mod
