module Cells
! 
! Author: A. C. Foster
! Puts trees on a 30x30 grid 
!
! Methods
!   initialize_cells


use Constants
use Tree

implicit none

!-------------------------------------------------------------------------------
! The Cells module contains the definition of the CellsData type, which
! holds attributes pertinent to all members of the cell type, and associated
! procedures.
!-------------------------------------------------------------------------------

  type                                              :: CellsData
       type(TreeData)                               :: tree
       logical                                      :: filled
  end type CellsData

contains

!-------------------------------------------------------------------------------
! Methods
!-------------------------------------------------------------------------------

  subroutine initialize_cells(self)

     class(CellsData),      intent(inout) :: self
     
     self%filled = .false.

  end subroutine initialize_cells

end module Cells
  




