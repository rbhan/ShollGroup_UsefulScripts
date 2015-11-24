#Supporting Functions for the Main Code

import numpy as np
import math as m
from mpi4py import MPI
COMM=MPI.COMM_WORLD   

####################################################################
def oneprint(myid,string):
     if myid == 0:
          print(string)
####################################################################
def find_line(myid,file,searchExp):
    if myid == 0:
        f=open(file,'r')
        index=0
        for line in f:
            if searchExp in line:
                break
            index += 1
        f.close()
    else:
         index = None
    
    index=COMM.bcast(index,root=0)   
    return(index)
####################################################################
def anyALTERED(iterable,string):
    for element in iterable:
        if element==string:
            return True
    return False
####################################################################
def cov_radius_assignment(atype,radii,atom):
     for i in range(0,len(atype)):
         if atom==atype[i,0]:
             cov_radius=radii[i,0]
             break
         elif i==len(atype):
             print('Unknown element:',atom)
     return(cov_radius)
####################################################################
def database_read(databases):
    file_covradii=open(databases+'/covalent_radii.txt','r')
    lines=file_covradii.readlines()
    file_covradii.close()
    
    atype=[]
    radii=[]
    for line in lines:
        elements=line.split()
        atype.append(elements[0])
        radii.append(elements[1])
        
    atype=np.array(atype,dtype='str')
    radii=np.array(radii,dtype='float')
    atype.shape=(len(atype),1)
    radii.shape=(len(radii),1)
    return(atype,radii)
####################################################################
def xyz_read(file_name):
    
    file_xyz=open(file_name,'r')
    lines=file_xyz.readlines()
    file_xyz.close() 
    
    #Store all the atom types and the xyz coordinates
    counter=0
    atom_type=[]
    xposi=[]
    yposi=[]
    zposi=[]
    for line in lines:
        counter=counter+1
        if counter > (2):
            column=line.split()
            atom_type.append(column[0])
            xposi.append(column[1])
            yposi.append(column[2])
            zposi.append(column[3])
            
    return(atom_type,xposi,yposi,zposi)
####################################################################
def transformation_matrix(aLEN,bLEN,cLEN,alphaD,betaD,gammaD,x_expand,y_expand,z_expand):
#Calcuating the transformation matrix 
#Fractional coordinates (a-axis is collinear with the x-axis)
    alpha_rad=alphaD*(np.pi/180)
    beta_rad=betaD*(np.pi/180)
    gamma_rad=gammaD*(np.pi/180)
    a=aLEN*x_expand
    b=bLEN*y_expand
    c=cLEN*z_expand
    
    c1 = c*m.cos(beta_rad)
    c2 = c*(m.cos(alpha_rad) - m.cos(gamma_rad) * m.cos(beta_rad))/m.sin(gamma_rad)
    c3 = (c**2 - c1**2 - c2**2)**0.5
    #C = [ a b*cos(gamma) c1 ; 0 b*sin(gamma) c2; 0 0 c3];
    C_mat=np.zeros([3,3])
    C_mat[0][0] = a
    C_mat[0][1] = b*m.cos(gamma_rad)
    C_mat[0][2] = c1
    C_mat[1][1] = b*m.sin(gamma_rad)
    C_mat[1][2] = c2
    C_mat[2][2] = c3
    C_mat_inv=np.linalg.inv(C_mat)
    #print C_mat
    #print C_mat_inv
    return(C_mat,C_mat_inv)
####################################################################
def connectivity_mapping(natoms,radii,atom_type,atype,coord_frac,C_mat,x_expand,y_expand,z_expand):
    rvec=np.zeros([3,1])
    rx=np.zeros([natoms,natoms])
    ry=np.zeros([natoms,natoms])
    rz=np.zeros([natoms,natoms])
    dummy=np.zeros([3,1])
    conn_matrix=np.zeros([natoms,13])
    conn_matrix_char=np.zeros([natoms,6],dtype=('a3 a3 a3 a3 a3 a3'))
    conn_matrix_nums=np.zeros([natoms,6],dtype=('a5 a5 a5 a5 a5 a5'))
    bond_fix = 0.4
    for i in range(natoms):
         index1=1
         index2=0
         atom1=atom_type[i]
         #This is extremely slow.... rad1=cov_radius_assignment(atom1)
         #c1=np.where(unique_list_sorted==atom1)
         #rad1=unique_rad_array[c1][0]
         rad1=cov_radius_assignment(atype,radii,atom1)
         conn_matrix[i][0]=i+1
         #print i
         #print i+1,rad1
         for j in range(natoms):
               atom2=atom_type[j]
               #This is extremely slow.... rad2=cov_radius_assignment(atom2)
               #c2=np.where(unique_list_sorted==atom2)
               #rad2=unique_rad_array[c2][0]
               rad2=cov_radius_assignment(atype,radii,atom2)
               #print j+1,rad2
           
               rx[i][j] = coord_frac[i,0] - coord_frac[j,0]
               ry[i][j] = coord_frac[i,1] - coord_frac[j,1]
               rz[i][j] = coord_frac[i,2] - coord_frac[j,2]
	
               #This includes PBCS...
               rx[i][j] = rx[i,j] - x_expand*round(rx[i,j]/float(x_expand))
               ry[i][j] = ry[i,j] - y_expand*round(ry[i,j]/float(y_expand))
               rz[i][j] = rz[i,j] - z_expand*round(rz[i,j]/float(z_expand))
	
               rvec[0][0] = rx[i,j]
               rvec[1][0] = ry[i,j]
               rvec[2][0] = rz[i,j]
               #print rvec
               rdummy = np.dot(C_mat,rvec) #Converting back to real coordinates
               #print rdummy
               r = (rdummy[0,0]**2 +  rdummy[1,0]**2 +  rdummy[2,0]**2)**0.5
               #print r 
               if i != j:
                    #print r,(rad1+rad2)
                    if r < 10e-8:
                         print('WARNING THERE ARE DUPLICATE ATOMS. REVISE XYZ FILE ACCORDINGLY.')
                         #print j+1
                         #print r
                    if r <= (rad1+rad2)+bond_fix: #This means atom(i) is bonded to atom(j)
                         #print j+1
                         #print atom1,rad1,atom2,rad2,r[i,j] 
                         conn_matrix[i][index1]=j+1
                         index1=index1+1
                         conn_matrix[i][index1]=r
                         index1=index1+1
                         conn_matrix_char[i][index2]=atom_type[j]
                         conn_matrix_nums[i][index2]=str(j+1)
                         #conn_matrix_convert[i][index2]=str(j+1)
                         index2=index2+1
                         #stag_array=np.array([coord_sorted[i][0],coord_sorted[i][1],coord_sorted[i][2],coord_sorted[j][0],coord_sorted[j][1],coord_sorted[j][2]])
                         #stag_array.shape=(1,6)
                         #bonds_4_c4d=np.concatenate((bonds_4_c4d,stag_array),axis=0)

    return(conn_matrix,conn_matrix_nums,conn_matrix_char)
####################################################################
def bond_maker(natoms,conn_matrix_nums):
    bond_array=np.empty([1,2],dtype='int')

    for i in range(natoms):
        for j in range(6):
            if i == 0 and j == 0 and conn_matrix_nums[i,j] != '': #start it out
                bond_array[0,0]=int(i+1)
                bond_array[0,1]=int(conn_matrix_nums[i,j])
            elif conn_matrix_nums[i,j] != '': #append the rest
                bond_array=np.append(bond_array,[[int(i+1),int(conn_matrix_nums[i,j])]],axis=0)

    bond_array=np.sort(bond_array)
    bond_array=unique_rows(bond_array)
    return(bond_array)
####################################################################
def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
####################################################################
def expand_LINKERPOS(natomsN,coord_fracN,x_expandN,y_expandN,z_expandN):

     #Then expand the actual unit cell
     xpos_newff1=[]
     ypos_newff1=[]
     zpos_newff1=[]

     linker_description_new=[]
     for ucxn in range(-1*x_expandN+1,x_expandN):
         for ucyn in range(-1*y_expandN+1,y_expandN):
             for uczn in range(-1*z_expandN+1,z_expandN):
                 for atomn in range(natomsN):
                      xpos_newff1.append(float(coord_fracN[atomn,0])+ucxn)
                      ypos_newff1.append(float(coord_fracN[atomn,1])+ucyn)
                      zpos_newff1.append(float(coord_fracN[atomn,2])+uczn)
                      linker_description_new.append(int(coord_fracN[atomn,3]))

     xpos_newff=xpos_newff1
     ypos_newff=ypos_newff1
     zpos_newff=zpos_newff1

     xpos_newff=np.array(xpos_newff,dtype='f')
     ypos_newff=np.array(ypos_newff,dtype='f')
     zpos_newff=np.array(zpos_newff,dtype='f')
     coord_fracnn=np.vstack([xpos_newff,ypos_newff,zpos_newff])

     linker_description_new=np.array(linker_description_new,dtype='int')

     return(coord_fracnn,linker_description_new)
