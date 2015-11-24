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
def transformation_matrix(xpos,ypos,zpos,aLEN,bLEN,cLEN,alphaD,betaD,gammaD):
#Calcuating the transformation matrix 
#Fractional coordinates (a-axis is collinear with the x-axis)
    alpha_rad=alphaD*(np.pi/180)
    beta_rad=betaD*(np.pi/180)
    gamma_rad=gammaD*(np.pi/180)
    #xpos_f=np.array(xpos,dtype='f')
    #ypos_f=np.array(ypos,dtype='f')
    #zpos_f=np.array(zpos,dtype='f')
    #a=np.max(np.abs(xpos_f))
    #b=np.max(np.abs(ypos_f))
    #c=np.max(np.abs(zpos_f))
    #all=np.array([a,b,c])
    #a=np.max(all)
    #b=np.max(all)
    #c=np.max(all)
    a=aLEN
    b=bLEN
    c=cLEN
    
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
############################################################################################################
def expand_UC(myid,natoms,atom_type,coord_frac,x_expand,y_expand,z_expand):
     #First expand the atom_type vector
     nunitcells=x_expand*y_expand*z_expand
     atom_type_new=[]
     for n in range(nunitcells):
         for atom in range(len(atom_type)):
             atom_type_new.append(atom_type[atom])

     #Then expand the actual unit cell
     lx=aLEN
     xytilt=bLEN*np.cos(np.radians(gammaD))
     xztilt=cLEN*np.cos(np.radians(betaD))
     ly=np.sqrt(bLEN**2-(xytilt**2))
     yztilt=(bLEN*cLEN*np.cos(np.radians(alphaD))-xytilt*xztilt)/float(ly)
     lz=np.sqrt(cLEN**2 - xztilt**2 - yztilt**2)
     print(lx,ly,lz)
     print(xytilt,xztilt,yztilt)

     xpos_newf=[]
     ypos_newf=[]
     zpos_newf=[]
     for ucx in range(x_expand):
         for ucy in range(y_expand):
             for ucz in range(z_expand):
                 for atom in range(natoms):
                      xpos_newf.append(float(coord_frac[atom,0])+ucx)
                      ypos_newf.append(float(coord_frac[atom,1])+ucy)
                      zpos_newf.append(float(coord_frac[atom,2])+ucz)

     xpos_newf=np.array(xpos_newf,dtype='f')
     ypos_newf=np.array(ypos_newf,dtype='f')
     zpos_newf=np.array(zpos_newf,dtype='f')
     coord_fracn=np.vstack([xpos_newf,ypos_newf,zpos_newf])

     return(atom_type_new,coord_fracn)
############################################################################################################ 

###Read xyz file
     oneprint(myid,'XYZ read.')
     atom_type,xpos,ypos,zpos=xyz_read('pim1.xyz') #do file format checking to catch errors
     natoms=len(atom_type)

  ###Expand the unit cell
     xposnpc=np.array(xpos,dtype='f')
     yposnpc=np.array(ypos,dtype='f')
     zposnpc=np.array(zpos,dtype='f')
     C_mat,C_mat_inv=transformation_matrix(xpos,ypos,zpos,12,12,25,90,90,90)
     coord_reg=np.vstack([xposnpc,yposnpc,zposnpc])
     coord_reg=coord_reg.T
     print(coord_reg.shape)
     coord_frac=np.dot(C_mat_inv,coord_reg.T)
     coord_frac=np.array(coord_frac.T)
     print(coord_frac.shape)
     #print(coord_frac)

     x_expand=1#-------------------------------------------------------------------------------------------------------------------> GLOBAL VARIABLE
     y_expand=1#-------------------------------------------------------------------------------------------------------------------> GLOBAL VARIABLE
     z_expand=1#-------------------------------------------------------------------------------------------------------------------> GLOBAL VARIABLE
     atom_type,coord_frac=expand_UC(myid,natoms,atom_type,coord_frac,x_expand,y_expand,z_expand)
     coord_frac=np.array(coord_frac.T)
     coord_reg=np.dot(C_mat,coord_frac.T)
     atom_type=np.array(atom_type,dtype='S7')
     natoms=len(atom_type)
     coord_reg=coord_reg.T
     print(coord_reg.shape)
     xpos=coord_reg[:,0]
     ypos=coord_reg[:,1]
     zpos=coord_reg[:,2]

     file=open('TEST_EXPAND.xyz','w')
     file.write(str(len(atom_type))+'\n')
     file.write(' '+'\n')
     for i in range(len(atom_type)):
         file.write(str(atom_type[i]).replace("['",'').replace("']",'')+' '+str(float(xpos[i]))+' '+str(float(ypos[i]))+' '+str(float(zpos[i]))+'\n')
     file.write('\r')
     file.close()