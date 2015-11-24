! file: helloworld.f90
! want to make this fortan code "inherently parallel"
! COMPILE LINE:::::f2py -c -m HELLLLO hellowworld.f90 --fcompiler=intelem -lmpi --fcompiler=intelem
!Currently Loaded Modulefiles:
  !1) python/2.7(default)   2) intel/13.2.146        3) openmpi/1.5.4         4) hwloc/1.2(default)    5) mkl/11.0

subroutine sayhello(comm,x)

  implicit none
  include '/usr/local/packages/openmpi/1.5.4/intel-13.2.146/include/mpif.h'
  integer, intent(inout) :: x
  !f2py intent(in,out) :: x

  integer :: comm !, rank, size, ierr
  !call MPI_Comm_size(comm, size, ierr) neither of these lines work, I want them to...
  !call MPI_Comm_rank(comm, rank, ierr) I really want to see which processor this is running on (is it on all or just one???)

  !do some stuff
  x=x+1 

  !print *,x
  !print *, 'Hello, World!'

end subroutine sayhello