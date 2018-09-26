
function cmp_function(a, b)
    use cpl_interface, only: segment
    type(segment), intent(in) :: a, b
    integer(2) :: cmp_function

    if (a.global_offset > b.global_offset) then
        cmp_function = 1
    elseif (b.global_offset > a.global_offset) then
        cmp_function = -1
    else
        cmp_function = 0
    endif

endfunction cmp_function

module qsort
! Define an overload of the default QSORT signature
! with a signature using the segment type.

  use ifport

  implicit none

interface
    subroutine qsort_segments(array, len, isize, comp)
       use cpl_interface, only: segment

       type(segment), dimension(len) :: array
       integer :: len, isize
       integer(2), external :: comp
       !
       ! Hook the overload to the real thing but be careful
       ! to connect to the correct qsort: the Fortran one, not
       ! the C one!
       !
       ! We need to call the _Fortran_ qsort, not the _C_ one, or
       ! there will be errors from the 1-origin vs. 0-origin indexing
       ! and the row-major vs. column-major ordering.
       !
       !DIR$ ATTRIBUTES ALIAS:'_qsort' :: qsort_segments
    end subroutine qsort_segments
end interface

endmodule qsort
