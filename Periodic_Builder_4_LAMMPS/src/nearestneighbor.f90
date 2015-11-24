subroutine NNlistgrid(comm,nnarray,xposC,yposC,zposC,xpos,ypos,zpos,pos,neg,Ngrid,Dim,NAtom)
    !include '/usr/local/packages/openmpi/1.5.4/intel-13.2.146/include/mpif.h'
    implicit none
    include '/usr/local/packages/openmpi/1.5.4/intel-13.2.146/include/mpif.h'
    integer, intent(in) :: Dim,NAtom,Ngrid
    real(8), intent(in) :: pos,neg

    real(8), intent(inout), dimension(0:NAtom-1,0:Dim-1) :: nnarray
    !f2py intent(in,out) :: nnarray

    real(8), intent(in), dimension(0:NAtom-1) :: xposC
    !f2py intent(in) :: xposC
    real(8), intent(in), dimension(0:NAtom-1) :: yposC
    !f2py intent(in) :: yposC
    real(8), intent(in), dimension(0:NAtom-1) :: zposC
    !f2py intent(in) :: zposC

    real(8), intent(in), dimension(0:Ngrid-1) :: xpos
    !f2py intent(in) :: xpos
    real(8), intent(in), dimension(0:Ngrid-1) :: ypos
    !f2py intent(in) :: ypos
    real(8), intent(in), dimension(0:Ngrid-1) :: zpos
    !f2py intent(in) :: zpos

    integer :: i,j,k,icounter
    real(8) :: dSQ

    integer :: comm 

    print *,Dim
    do i = 0, NAtom-1
        icounter = 0
        do j=0, Ngrid-1
            dSQ=(xposC(i)-xpos(j))**2+(yposC(i)-ypos(j))**2+(zposC(i)-zpos(j))**2

            if (dSQ .LE. pos .AND. dSQ .LE. neg) then
                icounter = icounter + 1
                nnarray(i,icounter-1)=j
                if (icounter == 6) then
                    exit 
                end if
            end if

            if (icounter /= 6) then
                do k = icounter, 6
                    nnarray(i,k)= -888888
                end do
            end if
        end do
    end do
    print *,'I AM DONE!!!'
end subroutine