subroutine bond_fix(BOND_ARRAY_sortedF,flag_arrayF,Nbonds,Dim2)
    implicit none

    integer, intent(in) :: Nbonds,Dim2

    integer, intent(in), dimension(0:Nbonds-1,0:Dim2-1) :: BOND_ARRAY_sortedF
    !f2py intent(in) :: BOND_ARRAY_sortedF

    integer, intent(inout), dimension(0:Nbonds-1,0) :: flag_arrayF
    !f2py intent(in,out) :: flag_arrayF

    integer :: i,j

    print *,'started'
    do i=0,nbonds
        do j=(i+1),nbonds
            if (BOND_ARRAY_sortedF(i,0) .eq. BOND_ARRAY_sortedF(j,0)) then
                if (BOND_ARRAY_sortedF(i,1) .eq. BOND_ARRAY_sortedF(j,1)) then
                    flag_arrayF(i,0)=1
                end if
            end if
        end do 
    end do
    print *,'ended'

end subroutine