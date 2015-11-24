subroutine connectivity_mapping(Natoms,atom_typeN,radii_arrayN,coord_fracN,C_matN,x_expandN,y_expandN,z_expandN,Conn_matrixN,Conn_matrix_numsN,bond_fixN,dima,dimb,dimc)
    implicit none

    integer, intent(in) :: Natoms,dima,dimb,dimc

    integer, intent(in) :: x_expandN,y_expandN,z_expandN
    !f2py intent(in) :: x_expandN,y_expandN,z_expandN
    real(8), intent(in) :: bond_fixN
    !f2py intent(in) :: bond_fixN

    integer, intent(in), dimension(0:Natoms-1) :: atom_typeN
    !f2py intent(in) :: atom_typeN
    real(8), intent(in), dimension(0:Natoms-1) :: radii_arrayN
    !f2py intent(in) :: radii_arrayN
    real(8), intent(in), dimension(0:Natoms-1,0:dima-1) :: coord_fracN
    !f2py intent(in) :: coord_fracN
    real(8), intent(in), dimension(0:dima-1,0:dima-1) :: C_matN
    !f2py intent(in) :: C_matN

    real(8), dimension(0:dima-1) :: rvecN
    real(8), dimension(0:dima-1) :: rdummyN
    real(8) :: rxN
    real(8) :: ryN
    real(8) :: rzN

    real(8), intent(inout), dimension(0:Natoms-1,0:dimc-1) :: Conn_matrixN
    !f2py intent(in,out) :: Conn_matrixN
    integer(8), intent(inout), dimension(0:Natoms-1,0:dimb-1) :: Conn_matrix_numsN
    !f2py intent(in,out) :: Conn_matrix_numsN

    integer :: i,j,index1,index2
    real(8) :: r,rad1,rad2
    character(8) :: atom1,atom2
    print *,'FORTRAN STARTED.'
    do i=0,natoms-1
         !print *,atom_typeN(i)
         !if (atom_typeN(i) .NE. 5) then
         !      print *,'test good'
         !endif   

         !print *,i
         index1=1
         index2=0
         rad1=radii_arrayN(i)
         Conn_matrixN(i,0)=i+1
         do j=0,Natoms-1
               rad2=radii_arrayN(j)
           
               rxN = coord_fracN(i,0) - coord_fracN(j,0)
               ryN = coord_fracN(i,1) - coord_fracN(j,1)
               rzN = coord_fracN(i,2) - coord_fracN(j,2)
	
               rxN = rxN - x_expandN*NINT(rxN/float(x_expandN))
               ryN = ryN - y_expandN*NINT(ryN/float(y_expandN))
               rzN = rzN - z_expandN*NINT(rzN/float(z_expandN))
	
               rvecN(0) = rxN
               rvecN(1) = ryN
               rvecN(2) = rzN
               rdummyN = MATMUL(C_matN,rvecN) 
               !print *,rdummy
               r = (rdummyN(0)**2 +  rdummyN(1)**2 +  rdummyN(2)**2)**0.5

               if (i /= j) then
                    if (r < 10e-8) then
                         print *,'WARNING THERE ARE DUPLICATE ATOMS. REVISE XYZ FILE ACCORDINGLY.'
                         print *,j+1
                         print *,r
                    end if
                    if (r <= rad1+rad2+bond_fixN) then 
                         !MAKE ZN ONLY BOND TO N
                         if (atom_typeN(i) .EQ. 5) then
                              if (atom_typeN(j) .EQ. 3) then
                                   Conn_matrixN(i,index1)=j+1
                                   index1=index1+1
                                   Conn_matrixN(i,index1)=r
                                   index1=index1+1
                                   Conn_matrix_numsN(i,index2)=j+1
                                   index2=index2+1
                              end if
                         !MAKE C NOT BOND TO ZN
                         else if (atom_typeN(i) .EQ. 2) then
                              if (atom_typeN(j) .NE. 5) then
                                   Conn_matrixN(i,index1)=j+1
                                   index1=index1+1
                                   Conn_matrixN(i,index1)=r
                                   index1=index1+1
                                   Conn_matrix_numsN(i,index2)=j+1
                                   index2=index2+1
                              end if
                         !MAKE N NOT BOND TO N
                         else if (atom_typeN(i) .EQ. 3) then
                              if (atom_typeN(j) .NE. 3) then
                                   Conn_matrixN(i,index1)=j+1
                                   index1=index1+1
                                   Conn_matrixN(i,index1)=r
                                   index1=index1+1
                                   Conn_matrix_numsN(i,index2)=j+1
                                   index2=index2+1
                              end if
                         !CATCH ALL
                         else
                              Conn_matrixN(i,index1)=j+1
                              index1=index1+1
                              Conn_matrixN(i,index1)=r
                              index1=index1+1
                              Conn_matrix_numsN(i,index2)=j+1
                              index2=index2+1
                         end if

                    end if
              end if
         end do
    end do
    print *,'FORTAN DONE.'
end subroutine

subroutine bond_fix(BOND_ARRAY_sortedF,flag_arrayF,Nbonds,Dim2)
    implicit none

    integer, intent(in) :: Nbonds,Dim2

    integer, intent(in), dimension(0:Nbonds-1,0:Dim2-1) :: BOND_ARRAY_sortedF
    !f2py intent(in) :: BOND_ARRAY_sortedF

    integer, intent(inout), dimension(0:Nbonds-1) :: flag_arrayF
    !f2py intent(in,out) :: flag_arrayF

    integer :: i,j

    print *,'started'
    do i=0,Nbonds
        do j=(i+1),Nbonds
            if (BOND_ARRAY_sortedF(i,0) .eq. BOND_ARRAY_sortedF(j,0)) then
                if (BOND_ARRAY_sortedF(i,1) .eq. BOND_ARRAY_sortedF(j,1)) then
                    flag_arrayF(i)=1
                end if
            end if
        end do 
    end do
    print *,'ended'

end subroutine