subroutine add_one(a,b,N) bind(C)
    implicit none
    integer(8),intent(in) :: N
    real(8),dimension(N),intent(in) :: a
    real(8),dimension(N),intent(inout):: b
    b(:) = a(:) + 1
   end subroutine add_one