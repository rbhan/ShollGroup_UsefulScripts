#Sholl Group / Phillips 66
#Purpose: LAMMPS Input,Connectivity Mapping    
#Name: Ross Verploegh
#GTaccountID: rverploegh3
#Date Originated: December 13, 2013
#Date Revised: December 13, 2013
############################################################################################################
#from mpi4py import MPI
#COMM=MPI.COMM_WORLD
import time
import os
import sys
import numpy as np
import itertools
import copy
from string import digits
from collections import defaultdict
import collections
import math

from support_functions import *
from sssr import *
import CONNECT
import BONDFIX
import COMBINEDv2
############################################################################################################
def read_inputs():
     if myid == 0:
          file=open(input+'/parameters.in','r')
          lines=file.readlines()
          file.close()

          #GLOBAL VARIABLES
          index=find_line(myid,input+'/parameters.in','SEED NUMBER')
          global GLOBAL_SEED
          GLOBAL_SEED = int(lines[index+1])

          index=find_line(myid,input+'/parameters.in','MAIN STRUCTURE FILENAME')
          global file_argument
          file_argument = str(lines[index+1]).replace('\n','')

          index=find_line(myid,input+'/parameters.in','EXPAND THE UNIT CELL')
          global x_expand
          global y_expand
          global z_expand
          elements=lines[index+1].split()
          x_expand= int(elements[0])
          y_expand= int(elements[1])
          z_expand= int(elements[2])

          index=find_line(myid,input+'/parameters.in','LJ CUTOFF')
          global lj_cutoff
          lj_cutoff = float(lines[index+1])

          index=find_line(myid,input+'/parameters.in','BOND CORRECTION FACTOR')
          global bond_corr
          bond_corr = float(lines[index+1])

          index=find_line(myid,input+'/parameters.in','METAL TO NEGLECT')
          global neglect_atype
          neglect_atype = str(lines[index+1]).replace('\n','')

          index=find_line(myid,input+'/parameters.in','NUMBER OF LINKER NEAREST NEIGHBORS')
          global n_linkerNNs
          n_linkerNNs = int(lines[index+1])

          index=find_line(myid,input+'/parameters.in','UNIT CELL LENGTHS')
          global aLEN
          global bLEN
          global cLEN
          global alphaD
          global betaD
          global gammaD
          elements=lines[index+1].split()
          aLEN = float(elements[0])
          bLEN = float(elements[1])
          cLEN = float(elements[2])
          elements=lines[index+2].split()
          alphaD = float(elements[0])
          betaD = float(elements[1])
          gammaD = float(elements[2])

          index=find_line(myid,input+'/parameters.in','TRICLINIC UNIT CELL?')
          global triclinic_flag
          elements=lines[index+1].split()
          triclinic_flag=str(elements[0])
          if alphaD != 90.0 or betaD!=90.0 or gammaD!=90.0: 
            if triclinic_flag!='YES':
              print('WARNING:: NOT POSSIBLE TO NOT HAVE A TRICLINIC BOX')
              time.sleep(100000)

          index=find_line(myid,input+'/parameters.in','FRAGMENT NAME')
          global partname
          partname = str(lines[index+1]).replace('\n','')

          index=find_line(myid,input+'/parameters.in','NUMBER OF ATOM IN OLD AND NEW LINKER')
          global natomsinlinker
          natomsinlinker=[]
          elements=lines[index+1].split()
          natomsinlinker.append(int(elements[0]))
          natomsinlinker.append(int(elements[1]))

          index=find_line(myid,input+'/parameters.in','PRIMARY CARBON AND NITROGEN LABELS FOR OLD LINKER')
          global primary_carbonORG
          global primary_nitrogenORG1
          global primary_nitrogenORG2
          elements=lines[index+1].split()
          primary_carbonORG=str(elements[0]).replace('\n','')
          primary_nitrogenORG1=str(elements[1]).replace('\n','')
          primary_nitrogenORG2=str(elements[2]).replace('\n','')

          index=find_line(myid,input+'/parameters.in','FRACTION OF LINKERS TO SWITCH')
          global frac2switch
          frac2switch= float(lines[index+1])

          index=find_line(myid,input+'/parameters.in','SHORT RANGE ORDER PARAMETER')
          global set_SOP
          global MCbeta
          global RMC_errorTOL
          elements=lines[index+1].split()
          set_SOP=float(elements[0])
          MCbeta=float(elements[1])
          RMC_errorTOL=float(elements[2])

          index=find_line(myid,input+'/parameters.in','INTENSITY CURVE INPUTS')
          global intensity_flag
          global MAJORITY_HYDROGEN
          global MINORITY_HYDROGEN
          elements=lines[index+1].split()
          intensity_flag=str(elements[0]).replace('\n','')
          MAJORITY_HYDROGEN=str(elements[1]).replace('\n','')
          MINORITY_HYDROGEN=str(elements[2]).replace('\n','')

          global DISTANCE_SPACING_TYPE
          elements=lines[index+2].split()
          DISTANCE_SPACING_TYPE=str(elements[0]).replace('\n','')

          global distances_to_look_fromD
          global KRISHNA_DATA
          if DISTANCE_SPACING_TYPE=='RANGE':
            global startD
            global stopD
            global stepD
            startD=int(elements[1])
            stopD =int(elements[2])
            stepD =float(elements[3])
            distances_to_look_fromD=np.arange(start=startD,stop=stopD,step=stepD)
          elif DISTANCE_SPACING_TYPE=='USER_DEFINED':
            global file_name4distances
            file_name4distances=str(elements[1]).replace('\n','')
            print(file_name4distances)
            distance_array=[]
            file=open(input+'/'+str(file_name4distances),'r')
            linestemp=file.readlines()
            for line in linestemp:
              ds=line.split()
              distance_array.append(float(ds[0]))
            file.close()
            distances_to_look_fromD=np.array(distance_array,dtype='float')
            print(distances_to_look_fromD)
           
            global file_name4data
            file_name4data=str(elements[2]).replace('\n','')
            print(file_name4data)
            data_array=[]
            file=open(input+'/'+str(file_name4data),'r')
            linestemp=file.readlines()
            for line in linestemp:
              ds=line.split()
              data_array.append(float(ds[1]))
            file.close()
            KRISHNA_DATA=np.array(data_array,dtype='float')
            print(KRISHNA_DATA)

          index=find_line(myid,input+'/parameters.in','NUMBER OF ATOM TYPES, BONDS, ANGLES, PROPER TORSIONS, IMPROPER TORSIONS')
          global natypes
          global nbondtypes
          global nangletypes
          global ntorsiontypes
          global nimpropertypes
          elements=lines[index+1].split()
          natypes = int(elements[0])
          nbondtypes = int(elements[1])
          nangletypes = int(elements[2])
          ntorsiontypes = int(elements[3])
          nimpropertypes= int(elements[4])

          index=find_line(myid,input+'/parameters.in','ATOM TYPES, CHARGES [e], LJ PARAMETERS [sigma A, epsilon kcal/mol]')
          global atypelist 
          global q
          global sigma 
          global epsilon
          global connection_def
          global atypedict
          global atypedictB
          global masses
          global atom_type_levels
          atypelist = []
          atypedict = {}
          atypedictB = {}
          q = []
          sigma = []
          epsilon =[]
          connection_def=[]
          masses=[]
          atom_type_levels=[]
          count=0
          for l in range(index+1,natypes+index+1):
              elements=lines[l].split()
              n=len(elements)
              atypelist.append(elements[0])
              q.append(elements[1])
              sigma.append(elements[2])
              epsilon.append(elements[3])
              masses.append(elements[4])
              atom_type_levels.append(elements[5])
              connection_def.append(elements[6:n])
              count+=1
              tempdict={elements[0]:count}
              atypedict.update(tempdict)
              tempdictB={count:elements[0]}
              atypedictB.update(tempdictB)

          print(atom_type_levels)
          #time.sleep(10000)

          index=find_line(myid,input+'/parameters.in','BOND K CONSTANT CORRECTION')
          global bond_KCORRECTION
          bond_KCORRECTION= float(lines[index+1])
          
          index=find_line(myid,input+'/parameters.in','BONDS (STRETCHING)')
          global bonddict
          global bond_param
          bonddict=defaultdict(list)
          bond_param = np.zeros((nbondtypes,2),dtype='float')
          for l in range(index+1,nbondtypes+index+1):
              elements=lines[l].split()
              bonddict[l-(index+1)+1]=[str(elements[0]),str(elements[1])]
              bond_param [l-(index+1),0]=float(elements[2])
              bond_param [l-(index+1),1]=float(elements[3])
          print(bonddict)

          index=find_line(myid,input+'/parameters.in','ANGLE K CONSTANT CORRECTION')
          global angle_KCORRECTION
          angle_KCORRECTION= float(lines[index+1])

          index=find_line(myid,input+'/parameters.in','ANGLES (BENDING)')
          global angledict
          global angle_param
          angledict=defaultdict(list)
          angle_param = np.zeros((nangletypes,2),dtype='float')
          for l in range(index+1,nangletypes+index+1):
              elements=lines[l].split()
              angledict[l-(index+1)+1]=[str(elements[0]),str(elements[1]),str(elements[2])]
              angle_param [l-(index+1),0]=float(elements[3])
              angle_param [l-(index+1),1]=float(elements[4])
          print(angledict)

          index=find_line(myid,input+'/parameters.in','PROPER TORSIONS [kcal/mol degree m]')
          global dihedraldict
          global torsion_param
          dihedraldict=defaultdict(list)
          torsion_param = np.zeros((ntorsiontypes,3),dtype='float')
          for l in range(index+1,ntorsiontypes+index+1):
              elements=lines[l].split()
              dihedraldict[l-(index+1)+1]=[str(elements[0]),str(elements[1]),str(elements[2]),str(elements[3])]
              torsion_param [l-(index+1),0]=float(elements[4])
              torsion_param [l-(index+1),1]=float(elements[5])
              torsion_param [l-(index+1),2]=float(elements[6])
          print(dihedraldict)

          index=find_line(myid,input+'/parameters.in','MODEL IMPROPERS AS PROPER')
          global improper_proper_flag
          elements=lines[index+1].split()
          improper_proper_flag=str(elements[0])

          index=find_line(myid,input+'/parameters.in','IMPROPER PLACEMENT POSITION')
          global centeratomINDEX
          elements=lines[index+1].split()
          centeratomINDEX=int(elements[0])

          index=find_line(myid,input+'/parameters.in','IMPROPER TORSIONS [kcal/mol degree m]')
          global improperdict
          global improper_param
          improperdict=defaultdict(list)
          improper_param = np.zeros((nimpropertypes,3),dtype='float')
          for l in range(index+1,nimpropertypes+index+1):
              elements=lines[l].split()
              improperdict[l-(index+1)+1]=[str(elements[0]),str(elements[1]),str(elements[2]),str(elements[3])]
              improper_param [l-(index+1),0]=float(elements[4])
              improper_param [l-(index+1),1]=float(elements[5])
              improper_param [l-(index+1),2]=float(elements[6])
          print(improperdict)

          index=find_line(myid,input+'/parameters.in','CHECK ATOM TYPING IN LAMMPS INPUT FILE?')
          global check_lmp_typing_flag
          elements=lines[index+1].split()
          check_lmp_typing_flag=str(elements[0])

     return(sigma,set_SOP)
############################################################################################################
def expand_UC(myid,natoms,atom_type,coord_frac,x_expand,y_expand,z_expand):
     #First expand the atom_type vector
     nunitcells=x_expand*y_expand*z_expand
     atom_type_new=[]
     for n in range(nunitcells):
         for atom in range(len(atom_type)):
             atom_type_new.append(atom_type[atom])

     #Then expand the actual unit cell
     global lx
     global ly
     global lz
     global xytilt
     global xztilt
     global yztilt
     lx=aLEN*x_expand
     xytilt=bLEN*y_expand*np.cos(np.radians(gammaD))
     xztilt=cLEN*z_expand*np.cos(np.radians(betaD))
     ly=np.sqrt((bLEN*y_expand)**2-(xytilt**2))
     yztilt=(bLEN*y_expand*cLEN*z_expand*np.cos(np.radians(alphaD))-xytilt*xztilt)/float(ly)
     lz=np.sqrt((cLEN*z_expand)**2 - xztilt**2 - yztilt**2)
     print(lx,ly,lz)
     print(xytilt,xztilt,yztilt)

     ALcheck=lx
     BLcheck=np.sqrt(ly**2 + xytilt**2)
     CLcheck=np.sqrt(lz**2+xztilt**2+yztilt**2)
     AAcheck=np.degrees(np.arccos(float(xytilt*xztilt+ly*yztilt)/float(ALcheck*BLcheck)))
     BAcheck=np.degrees(np.arccos(float(xztilt)/float(CLcheck)))
     CAcheck=np.degrees(np.arccos(float(xytilt)/float(BLcheck)))
     print(aLEN*x_expand,bLEN*y_expand,cLEN*z_expand,alphaD,betaD,gammaD)
     print(ALcheck,BLcheck,CLcheck,AAcheck,BAcheck,CAcheck)
     #time.sleep(10000)
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
def remove_all(element, list):
    return filter(lambda x: x != element, list)

def neighborscan(myid,NNlist,neighindex,NEGLECT_METAL_LIST):
    #if myid == 0:
    NNs=copy.deepcopy(NNlist[neighindex,:])
    #print(NNs)
    NNs=remove_all('',NNs)
    #I NEED TO REMOVE ALL THE METAL ATOM INDEXES AS WELL
    #print(NNs)
    kNNs=[]
    for n in range(len(NNs)):
        if (float(NNs[n])-1) < float(neighindex):
            kNNs.append(int(NNs[n])-1)
    #print(kNNs)
    kNNs=remove_metals(myid,NEGLECT_METAL_LIST,kNNs)
    #print(kNNs)
    kNNs=np.array(kNNs,dtype='int')
    #else:
        #kNNs=None

    #kNNs=MPI.COMM_WORLD.bcast(kNNs,root=0)                
    return(kNNs)

def metals_indexes(myid,natoms,ATOM_TYPE_REAL,neglect_atype):
    NEGLECT_METAL_LIST=[]
    for atom in range(natoms):
        if neglect_atype == ATOM_TYPE_REAL[atom]:
            NEGLECT_METAL_LIST.append(atom)
    return(NEGLECT_METAL_LIST)

def remove_metals(myid,NEGLECT_METAL_LIST,kNNs):
    NNs = [x for x in kNNs if x not in NEGLECT_METAL_LIST]
    return(NNs)
            
def findroot(myid,inti,ptr):
    #print('INSIDE')
    #print(inti)
    #print(ptr)
    while ptr[inti] > 0:
        inti = int(ptr[inti])
    return(inti)

def percolate(myid,neglect_atype,ATOM_TYPE_REAL,conn_matrix_nums,NEGLECT_METAL_LIST):  
    ptr=['E' for x in range(len(ATOM_TYPE_REAL))] 
    big=0
    #order=np.arange(len(sp_rad)) 
    #print('starting')
    for site in range(len(ATOM_TYPE_REAL)):
        #print(str(ATOM_TYPE_REAL[site]).replace('[','').replace(']',''))
        #print('SITE ::'+str(site))
        if str(ATOM_TYPE_REAL[site,0]) != str(neglect_atype):
            #print('blah')
            r1 = int(site)
            s1 = int(site)
            ptr[s1] = int(-1)
            #print('blah'+str(ptr[s1]))
            kNNs=neighborscan(myid,conn_matrix_nums,int(site),NEGLECT_METAL_LIST)
            #print('PASSED NN SCAN')
            #print(kNNs)
            for nNNs in range(len(kNNs)):
                s2 = int(kNNs[nNNs])
                #print('s2 '+str(s2))
                if str(ptr[s2]) != 'E': #this makes sure a NN is not empty
                    r2 = findroot(myid,s2,ptr)
                    #print('r2 '+str(r2))
                    if r2 != r1:
                        if int(ptr[r1])>int(ptr[r2]):
                            ptr[r2] += int(ptr[r1])
                            ptr[r1] = r2 
                            r1 = r2
                        else:
                            ptr[r1] += int(ptr[r2])
                            ptr[r2] = r1  

                        if -ptr[r1] > big:
                            big = int(-ptr[r1])   
        elif str(ATOM_TYPE_REAL[site,0]) == str(neglect_atype):
            ptr[site]=-1                       
                
    #print('big '+str(big))
    nclustersizes = [x for x in ptr if x < -1] 
    nclusters = len(nclustersizes) 
    nclusters=np.array(nclusters,dtype='int')
    NCL=copy.deepcopy(nclusters)
    NCL=np.array(NCL,dtype='int')
    NCL.resize(1,1)
    return(NCL,nclustersizes,ptr)

def linker_dict_creator(myid,natoms,nclustersizes,natomsinlinker,ptr):
    organic_placement=[]
    metal_placement=[]
    #print(ptr)
    #print(-natomsinlinker[0])
    for atom in range(natoms):
        if ptr[atom] != 'E':
            if (int(ptr[atom])) == -natomsinlinker[0]:
                organic_placement.append(atom)
            elif (int(ptr[atom])) == -natomsinlinker[1]:
                organic_placement.append(atom)

            if (int(ptr[atom])) == -1:
                metal_placement.append(atom)

    #print(organic_placement)
    #print(len(organic_placement))

    #print(metal_placement)
    #print(len(metal_placement))

    organic_dict={}
    for clust in range(len(nclustersizes)):
        tempdict={organic_placement[clust]:clust+1}
        organic_dict.update(tempdict)

    #od = collections.OrderedDict(sorted(organic_dict.items()))
    #print(od)
    #time.sleep(5)

    cluster_dict=defaultdict(list)
    for atom in range(natoms):  
        #print(atom)
        atom_root=findroot(myid,atom,ptr)
        #print(atom_root)
        if atom_root in organic_placement:
            #print('IN LIST')
            #if atom_root == 3588:
            #    print(CONNECT)
            clusterid=organic_dict[atom_root]
            #print('ASSIGNED')
            cluster_dict[clusterid].append(atom)
            #print('APPENDED')

    print('FINISHED!!!!!!')
    #print(cluster_dict)
    return(cluster_dict)
############################################################################################################
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return (vector / np.linalg.norm(vector))

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return (angle)

def midpoint(v1,v2):
    x=(v1[0]-v2[0])/float(2)
    y=(v1[1]-v2[1])/float(2)
    z=(v1[2]-v2[2])/float(2)
    vec=(x,y,z)
    return(vec)

def plane2normal(pt1,pt2,pt3):
    v1=[pt1[0]-pt2[0],pt1[1]-pt2[1],pt1[2]-pt2[2]]
    v2=[pt1[0]-pt3[0],pt1[1]-pt3[1],pt1[2]-pt3[2]]
    #print(v1)
    #print(v2)
    crossvec=np.cross(v1,v2)
    return(crossvec)

def axisfinder(v1,v2):
    crossvec=np.cross(v1,v2)
    return (crossvec/ np.linalg.norm(crossvec))

def YPR_transform(x,y,z,pt1,pt2,angle):
    acos=np.cos(np.radians(angle))
    asin=np.sin(np.radians(angle))
    d=pt1[0]
    e=pt1[1]
    f=pt1[2]
    u=pt2[0]
    v=pt2[1]
    w=pt2[2]
    unvec=unit_vector((d-u,e-v,f-w))
    a=unvec[0]
    b=unvec[1]
    c=unvec[2]

    #TRANSFORMATION MATRIX
    xpos=(a*(v**2+w**2) - u*(b*v + c*w - u*x - v*y - w*z))*(1-acos) + x*acos + (-c*v + b*w - w*y + v*z)*asin
    ypos=(b*(u**2+w**2) - v*(a*u + c*w - u*x - v*y - w*z))*(1-acos) + y*acos + (c*u - a*w + w*x - u*z)*asin
    zpos=(c*(u**2+v**2) - w*(a*u + b*v - u*x - v*y - w*z))*(1-acos) + z*acos + (-b*u + a*v - v*x + u*y)*asin

    return(xpos,ypos,zpos)

def YPR_transform2(x,y,z,pt1,pt2,angle):
    acos=np.cos(np.radians(angle))
    asin=np.sin(np.radians(angle))
    a=pt1[0]
    b=pt1[1]
    c=pt1[2]
    d=pt2[0]
    e=pt2[1]
    f=pt2[2]
    unvec=unit_vector((a-d,b-e,c-f))
    u=unvec[0]
    v=unvec[1]
    w=unvec[2]

    #TRANSFORMATION MATRIX
    xpos=(a*(v**2+w**2) - u*(b*v + c*w - u*x - v*y - w*z))*(1-acos) + x*acos + (-c*v + b*w - w*y + v*z)*asin
    ypos=(b*(u**2+w**2) - v*(a*u + c*w - u*x - v*y - w*z))*(1-acos) + y*acos + (c*u - a*w + w*x - u*z)*asin
    zpos=(c*(u**2+v**2) - w*(a*u + b*v - u*x - v*y - w*z))*(1-acos) + z*acos + (-b*u + a*v - v*x + u*y)*asin

    return(xpos,ypos,zpos)

def center_pos(xpos,ypos,zpos):
     xCENTER=np.sum(xpos)/float(len(xpos))
     yCENTER=np.sum(ypos)/float(len(ypos))
     zCENTER=np.sum(zpos)/float(len(zpos))

     return(xCENTER,yCENTER,zCENTER)
############################################################################################################
def types_charges(myid,natoms,atom_type,conn_matrix_nums,conn_matrix_char):
    NUMSTOAVOID=set('0123456789')
  ##LEVEL 1
    #This makes the dictionary
    permdict={}
    for x in range(natypes):
        if connection_def[x]=='NONE' or len(connection_def[x])==1:
            perm_list=''.join(connection_def[x])
        else:
            string_link=''.join(connection_def[x])
            perm_list=list(map("".join, itertools.permutations(string_link)))
        #print(perm_list)
        tempdict={(x+1):perm_list}
        permdict.update(tempdict)

    print(permdict[7])

    #conn_matrix_char_TYPED=np.zeros([natoms,6],dtype=('a3')) 
    conn_matrix_char_TYPED=copy.deepcopy(conn_matrix_char)
    LAMMPS_CHARGES=np.ones([natoms,1],dtype='f') 
    LAMMPS_ATOMINDEXES=np.ones([natoms,1],dtype='int') 
    ATOM_TYPE_REAL=np.ones([natoms,1],dtype='S7')
    NEGLECT_ATOM_INDEXES=[]  
    #print(atypelist)
    for atom in range(natoms):
        string_link=str(conn_matrix_char[atom][0])+str(conn_matrix_char[atom][1])+str(conn_matrix_char[atom][2])+\
                    str(conn_matrix_char[atom][3])+str(conn_matrix_char[atom][4])+str(conn_matrix_char[atom][5])        
        ATOMINDEX=1
        for i in range(natypes):
            if atom_type[atom] in atypelist:
                if permdict[atypedict[atom_type[atom]]]=='NONE':
                     ATOM_TYPE_REAL[atom,0]=atom_type[atom]
                     LAMMPS_CHARGES[atom,0]=q[atypedict[atom_type[atom]]-1]
                     LAMMPS_ATOMINDEXES[atom,0]=atypedict[atom_type[atom]]
                     break

            elif string_link in permdict[ATOMINDEX] and any((c in NUMSTOAVOID) for c in permdict[ATOMINDEX])==False:
                ATOM_TYPE_REAL[atom,0]=atypedictB[ATOMINDEX]
                LAMMPS_CHARGES[atom,0]=q[atypedict[atypedictB[ATOMINDEX]]-1]
                LAMMPS_ATOMINDEXES[atom,0]=atypedict[atypedictB[ATOMINDEX]]
                break
            ATOMINDEX+=1

  ##Fix according to LEVEL 1 screen
    #print(conn_matrix_nums)
    for atom in range(natoms):
        for x in range(6):
            if conn_matrix_nums[atom][x] != 0:
                aindex=int(conn_matrix_nums[atom][x])-1
                if ATOM_TYPE_REAL[aindex]!='1':
                    conn_matrix_char_TYPED[atom][x]=ATOM_TYPE_REAL[aindex,0]        
  
  ##LEVEL 2
    for atom in range(natoms):
        string_link=str(conn_matrix_char_TYPED[atom][0])+str(conn_matrix_char_TYPED[atom][1])+str(conn_matrix_char_TYPED[atom][2])+\
                    str(conn_matrix_char_TYPED[atom][3])+str(conn_matrix_char_TYPED[atom][4])+str(conn_matrix_char_TYPED[atom][5])
        #print(string_link)
        
        ATOMINDEX=1
        for i in range(natypes):
            if string_link in permdict[ATOMINDEX] and ATOM_TYPE_REAL[atom,0]=='1' and atom_type_levels[i]=='L2':
                #print(atom_type_levels[i])
                ATOM_TYPE_REAL[atom,0]=atypedictB[ATOMINDEX]
                LAMMPS_CHARGES[atom,0]=q[atypedict[atypedictB[ATOMINDEX]]-1]
                LAMMPS_ATOMINDEXES[atom,0]=atypedict[atypedictB[ATOMINDEX]]
                break
            ATOMINDEX+=1

  ##Fix according to LEVEL 2 screen
    for atom in range(natoms):
        for x in range(6):
            if conn_matrix_nums[atom][x] != 0:
                aindex=int(conn_matrix_nums[atom][x])-1
                if ATOM_TYPE_REAL[aindex]!='1':
                    conn_matrix_char_TYPED[atom][x]=ATOM_TYPE_REAL[aindex,0]

    #print(conn_matrix_char_TYPED)
    #time.sleep(10000)

  ##LEVEL 3
    for atom in range(natoms):
        string_link=str(conn_matrix_char_TYPED[atom][0])+str(conn_matrix_char_TYPED[atom][1])+str(conn_matrix_char_TYPED[atom][2])+\
                    str(conn_matrix_char_TYPED[atom][3])+str(conn_matrix_char_TYPED[atom][4])+str(conn_matrix_char_TYPED[atom][5])
        #print(string_link)
        
        ATOMINDEX=1
        for i in range(natypes):
            if string_link in permdict[ATOMINDEX] and ATOM_TYPE_REAL[atom,0]=='1' and atom_type_levels[i]=='L3':
                #print(atom_type_levels[i])
                ATOM_TYPE_REAL[atom,0]=atypedictB[ATOMINDEX]
                LAMMPS_CHARGES[atom,0]=q[atypedict[atypedictB[ATOMINDEX]]-1]
                LAMMPS_ATOMINDEXES[atom,0]=atypedict[atypedictB[ATOMINDEX]]
                break
            ATOMINDEX+=1

  ##Fix according to LEVEL 3 screen
    for atom in range(natoms):
        for x in range(6):
            if conn_matrix_nums[atom][x] != 0:
                aindex=int(conn_matrix_nums[atom][x])-1
                if ATOM_TYPE_REAL[aindex]!='1':
                    conn_matrix_char_TYPED[atom][x]=ATOM_TYPE_REAL[aindex,0]

  ##LEVEL 4
    for atom in range(natoms):
        string_link=str(conn_matrix_char_TYPED[atom][0])+str(conn_matrix_char_TYPED[atom][1])+str(conn_matrix_char_TYPED[atom][2])+\
                    str(conn_matrix_char_TYPED[atom][3])+str(conn_matrix_char_TYPED[atom][4])+str(conn_matrix_char_TYPED[atom][5])
        #print(string_link)
        
        ATOMINDEX=1
        for i in range(natypes):
            if string_link in permdict[ATOMINDEX] and ATOM_TYPE_REAL[atom,0]=='1' and atom_type_levels[i]=='L4':
                #print(atom_type_levels[i])
                ATOM_TYPE_REAL[atom,0]=atypedictB[ATOMINDEX]
                LAMMPS_CHARGES[atom,0]=q[atypedict[atypedictB[ATOMINDEX]]-1]
                LAMMPS_ATOMINDEXES[atom,0]=atypedict[atypedictB[ATOMINDEX]]
                break
            ATOMINDEX+=1

  ##Fix according to LEVEL 4 screen
    for atom in range(natoms):
        for x in range(6):
            if conn_matrix_nums[atom][x] != 0:
                aindex=int(conn_matrix_nums[atom][x])-1
                if ATOM_TYPE_REAL[aindex]!='1':
                    conn_matrix_char_TYPED[atom][x]=ATOM_TYPE_REAL[aindex,0]

    #print(ATOM_TYPE_REAL)
    #print(LAMMPS_CHARGES)
    #print(LAMMPS_ATOMINDEXES)
    #print(ATOM_TYPE_REAL[104])
    #print(int(conn_matrix_nums[275][2])-1)
    #print(conn_matrix_nums)         
    #print(conn_matrix_char[95][:])         
    #print(conn_matrix_char_TYPED[95][:])

    return(ATOM_TYPE_REAL,LAMMPS_ATOMINDEXES,LAMMPS_CHARGES,conn_matrix_char_TYPED)
############################################################################################################
def bonds_angles(myid,natoms,conn_matrix_char,conn_matrix_nums):
    print('part 1')
    BOND_A=[]
    BOND_B=[]
    ANGLE_A=[]
    ANGLE_B=[]
    ANGLE_C=[]
 
    bond_count=0
    angle_count=0   
    for atom in range(natoms):
        length_index=0
        for j in range(0,6):
            if conn_matrix_char[atom][j]!='':
                length_index+=1

        for k in range(length_index):
            bond_count+=1
            BOND_A.append(int(atom+1))
            BOND_B.append(int(conn_matrix_nums[atom][k]))
            if length_index >= 2:
                for l in range (k+1,length_index):
                    angle_count+=1
                    ANGLE_A.append(int(conn_matrix_nums[atom][k]))
                    ANGLE_B.append(int(atom+1))
                    ANGLE_C.append(int(conn_matrix_nums[atom][l]))

    bond_count=(bond_count/2) #This is need since there are bond duplicates
    BOND_A=np.array(BOND_A,dtype='int')
    BOND_B=np.array(BOND_B,dtype='int')
    ANGLE_A=np.array(ANGLE_A,dtype='int')
    ANGLE_B=np.array(ANGLE_B,dtype='int')
    ANGLE_C=np.array(ANGLE_C,dtype='int')

    BOND_ARRAY=np.vstack((BOND_A,BOND_B))
    ANGLE_ARRAY=np.vstack((ANGLE_A,ANGLE_B,ANGLE_C))

    BOND_ARRAY=BOND_ARRAY.transpose()
    ANGLE_ARRAY=ANGLE_ARRAY.transpose()

    #There are duplicates in the bond array so I need to remove them
    print('part 2')
    count=1
    BOND_ARRAY_sorted=copy.deepcopy(np.sort(BOND_ARRAY)) #sort the array
    flag_array=np.zeros([len(BOND_ARRAY)]) #set flags for bonds that need to be removed
    BOND_ARRAY_SHORT=np.zeros([len(BOND_ARRAY)/2,2])

    flag_array=COMBINEDv2.bond_fix(BOND_ARRAY_sorted,flag_array) #######THIS IS FUCKKKKKED################
#__________________________________________________
    #for i in range(len(BOND_ARRAY_sorted)):
        #for j in range(i+1,len(BOND_ARRAY_sorted)):
            #if BOND_ARRAY_sorted[i][0]==BOND_ARRAY_sorted[j][0]:
                 #if BOND_ARRAY_sorted[i][1]==BOND_ARRAY_sorted[j][1]:
                      #flag_array[i][0]=1
#__________________________________________________

    counterblah = -1
    for i in range(len(BOND_ARRAY_sorted)):
        if str(flag_array[i])=='0':
        #if flag_array[i][0]==0:
            counterblah+=1
            BOND_ARRAY_SHORT[counterblah][0]=int(BOND_ARRAY[i][0])
            BOND_ARRAY_SHORT[counterblah][1]=int(BOND_ARRAY[i][1])

    print('part 3')
    BOND_ARRAYback=BOND_ARRAY_SHORT
    #print(len(BOND_ARRAY))
    #print(len(BOND_ARRAY_SHORT))
    BOND_ARRAY_sortedback=BOND_ARRAY_SHORT
    BOND_ARRAY_sortedback=np.sort(BOND_ARRAY_sortedback)
    BOND_ARRAYback=BOND_ARRAYback.astype(int)
    BOND_ARRAYback=BOND_ARRAY_sortedback.astype(int)
    return(bond_count,angle_count,BOND_ARRAYback,ANGLE_ARRAY)
        
############################################################################################################
def dihedrals(myid,natoms,conn_matrix_char,conn_matrix_nums,BOND_ARRAY):
  dihedral_count=0
  DIHED_A=[]
  DIHED_B=[]
  DIHED_C=[]
  DIHED_D=[]

  for i in range(0,len(BOND_ARRAY)):
    lenx=0
    leny=0
    x=BOND_ARRAY[i][0]-1
    y=BOND_ARRAY[i][1]-1
    for j in range(0,4):
            if conn_matrix_char[x][j]!='':
                lenx=lenx+1
            if conn_matrix_char[y][j]!='':
                leny=leny+1
    for k in range(0,lenx):
        for l in range(0,leny):
            if conn_matrix_nums[x][k]!=(int(y+1)) and conn_matrix_nums[y][l]!=(int(x+1)):
                dihedral_count=dihedral_count+1
                DIHED_A.append(int(conn_matrix_nums[x][k]))
                DIHED_B.append(x+1)                
                DIHED_C.append(y+1)               
                DIHED_D.append(int(conn_matrix_nums[y][l]))

  DIHEDRAL_ARRAY=np.vstack((DIHED_A,DIHED_B,DIHED_C,DIHED_D))
  DIHEDRAL_ARRAY=DIHEDRAL_ARRAY.transpose()

  file=open('DIHEDRALS','w')
  for i in range(len(DIHEDRAL_ARRAY)):
    file.write(str(DIHEDRAL_ARRAY[i,0])+' '+str(DIHEDRAL_ARRAY[i,1])+' '+str(DIHEDRAL_ARRAY[i,2])+' '+str(DIHEDRAL_ARRAY[i,3])+'\n')
  file.close()
  
  #Remove duplicate dihedrals #I do not think this is even necessary
  flag_array=np.zeros([len(DIHEDRAL_ARRAY),1],dtype='int') #set flags for dihedrals that need to be removed
  #print(len(DIHEDRAL_ARRAY))
  #for i in range(len(DIHEDRAL_ARRAY)):
  #    print(i)
  #    temp1=np.sort([DIHEDRAL_ARRAY[i,1],DIHEDRAL_ARRAY[i,2]])
  #    temp1=str(temp1)
  #    for j in range(i+1,len(DIHEDRAL_ARRAY)):
  #      if i!=j:
  #        temp2=np.sort([DIHEDRAL_ARRAY[j,1],DIHEDRAL_ARRAY[j,2]])
  #        temp2=str(temp2)
  #        if cmp(temp1,temp2)==0:
  #            if   (DIHEDRAL_ARRAY[i,0]==DIHEDRAL_ARRAY[j,0]) and (DIHEDRAL_ARRAY[i,3]==DIHEDRAL_ARRAY[j,3]):
  #              flag_array[i,0]=1
  #            elif (DIHEDRAL_ARRAY[i,0]==DIHEDRAL_ARRAY[j,3]) and (DIHEDRAL_ARRAY[i,3]==DIHEDRAL_ARRAY[j,0]):
  #              flag_array[i,0]=1

  DIHEDRAL_ARRAY_SHORT=np.zeros([len(DIHEDRAL_ARRAY),4],dtype='int')
  counterblah = -1
  for i in range(len(DIHEDRAL_ARRAY)):
    if flag_array[i]==0:
      counterblah+=1
      DIHEDRAL_ARRAY_SHORT[counterblah][0]=int(DIHEDRAL_ARRAY[i][0])
      DIHEDRAL_ARRAY_SHORT[counterblah][1]=int(DIHEDRAL_ARRAY[i][1])
      DIHEDRAL_ARRAY_SHORT[counterblah][2]=int(DIHEDRAL_ARRAY[i][2])
      DIHEDRAL_ARRAY_SHORT[counterblah][3]=int(DIHEDRAL_ARRAY[i][3])
  dihedral_count=len(DIHEDRAL_ARRAY_SHORT)

  return(dihedral_count,DIHEDRAL_ARRAY_SHORT)
############################################################################################################
def ljparameters(epsilon,sigma):
    print('Lortentz-Berthelot Mixing Rules Calculator for LAMMPS')
    print('Check the input epsilon [K] and sigma [A] values.')

    #TRAPPE PARAMETERS FOR N-ALKANES
    #methane: 148.0 3.73 (CH4)
    #ethane: 98.0 3.75 (CH3) #both carbons are treated the same
    #propane: 98.0 3.750 (CH3)/46.0 3.950 (CH2)
    #butane: 98.0 3.750 (CH3)/46.0 3.950 (CH2)
    #isobutane: 98.0 3.750 (CH3)/ 10.0 4.680 (CH)
    #isobutene: 85.0 3.675 (CH2)/ 20.0 3.85 (C)/ 98.0 3.75 (CH3)
    #1-butene: 85.0 3.675 (CH2)/ 47.0 3.73 (CH)/ 46.0 3.950 (CH2)/ 98.0 3.750 (CH3) 
    #propene: 85.0 3.675 (CH2)/ 47.0 3.73 (CH)/ 98.0 3.75 (CH3)
    #ethene: 85.0 3.675

    file=open(output+'/VDW_LAMMPS.data','w')
    for i in range(natypes):
        for j in range(i,natypes):
            mixed_epsi=np.sqrt(float(epsilon[i])*float(epsilon[j]))#*0.00198717
            mixed_sigma=(float(sigma[i])+float(sigma[j]))/float(2)
            file.write('pair_coeff '+str(i+1)+' '+str(j+1)+' lj/cut '+str(mixed_epsi)+' '+str(mixed_sigma)+'\n')
    file.close()
    return()

############################################################################################################
def xyz_write(myid,name,natoms,ATOM_TYPE_REAL,coord_reg):
    file=open(output+'/'+str(name)+'TYPED.xyz','w')
    file.write(str(natoms)+'\n')
    file.write(str(name)+' '+str(aLEN)+' '+str(bLEN)+' '+str(cLEN)+' '+str(alphaD)+' '+str(betaD)+' '+str(gammaD)+'\n')
    for i in range(natoms):
        file.write(str(ATOM_TYPE_REAL[i]).replace("['",'').replace("']",'')+' '+str(coord_reg[i,0])+' '+str(coord_reg[i,1])+' '+str(coord_reg[i,2])+'\n')
    file.write('\r')
    file.close()

    file=open(output+'/'+str(name)+'UNTYPED.xyz','w')
    file.write(str(natoms)+'\n')
    file.write(str(name)+' '+str(aLEN)+' '+str(bLEN)+' '+str(cLEN)+' '+str(alphaD)+' '+str(betaD)+' '+str(gammaD)+'\n')
    for i in range(natoms):
        file.write(str(ATOM_TYPE_REAL[i]).translate(None,digits).replace("['",'').replace("']",'')+' '+str(coord_reg[i,0])+' '+str(coord_reg[i,1])+' '+str(coord_reg[i,2])+'\n')
    file.write('\r')
    file.close()
    print('Typed and Untyped XYZ files have been written.')
    return(1)
############################################################################################################
#BOND TYPING
def bondtyping(bonddict,BOND_ARRAY,atom_type,coord_regblah):
  fileBONDS=open(output+'/BONDS_TYPED_LENGTHS.data','w')

  REALBONDARRAY=np.ones([0,2],dtype='int')
  bondtypes=[]
  bond_lengths=[]
  bond_infoA=np.empty([len(bonddict),2],dtype='a3')
  bond_infoB=np.zeros([len(bonddict),4],dtype='float')

  for b in range(len(bonddict)):
    temp=bonddict[b+1]
    bstr1=temp[0:2]
    lmptype=int(b+1)
    for c in range(len(BOND_ARRAY)):
      bstr2=np.sort([atom_type[BOND_ARRAY[c,0]-1],atom_type[BOND_ARRAY[c,1]-1]])
      bstr2=[str(bstr2[0][0]),str(bstr2[1][0])]
      if set(bstr1)==set(bstr2):
        bondtypes.append(int(lmptype))
        bond_length=np.sqrt((coord_regblah[BOND_ARRAY[c,0]-1,0]-coord_regblah[BOND_ARRAY[c,1]-1,0])**2+(coord_regblah[BOND_ARRAY[c,0]-1,1]-coord_regblah[BOND_ARRAY[c,1]-1,1])**2+(coord_regblah[BOND_ARRAY[c,0]-1,2]-coord_regblah[BOND_ARRAY[c,1]-1,2])**2)
        if bond_length <= 10:
          fileBONDS.write(str(lmptype)+' '+str(bond_length)+'\n')
          bond_lengths.append(bond_length)
          list_bond_atom_types=bonddict[lmptype]
          bond_infoA[lmptype-1,0]=str(list_bond_atom_types[0])
          bond_infoA[lmptype-1,1]=str(list_bond_atom_types[1])
          bond_infoB[lmptype-1,0]=lmptype
          bond_infoB[lmptype-1,1]=bond_infoB[lmptype-1,1]+bond_length
          bond_infoB[lmptype-1,2]=bond_infoB[lmptype-1,2]+1
          bond_infoB[lmptype-1,3]=float(bond_infoB[lmptype-1,1])/float(bond_infoB[lmptype-1,2])

        REALBONDARRAY=np.vstack((REALBONDARRAY,[BOND_ARRAY[c,0],BOND_ARRAY[c,1]]))
  print('REAL NUMBER OF BONDS: '+str(len(bondtypes)))
  print('SIZE OF CORRECT BOND ARRAY: '+str(len(REALBONDARRAY)))

  fileBONDS.close()
  
  fileBONDS=open(output+'/SUMMARY_BONDS_TYPED_LENGTHS.data','w')
  for b in range(len(bonddict)):
    fileBONDS.write(str(bond_infoA[b,0])+' '+str(bond_infoA[b,1])+' '+str(int(bond_infoB[b,0]))+' '+str(bond_infoB[b,1])+' '+str(int(bond_infoB[b,2]))+' '+str(bond_infoB[b,3])+'\n')      
  fileBONDS.close()

  return(bondtypes,REALBONDARRAY)

#ANGLE TYPING
def angletyping(angledict,ANGLE_ARRAY,atom_type,coord_regblah):
  fileANGLES=open(output+'/ANGLES_TYPED_LENGTHS.data','w')

  REALANGLEARRAY=np.ones([0,3],dtype='int')
  angletypes=[]
  angle_infoA=np.zeros([len(angledict),3],dtype='a3')
  angle_infoB=np.zeros([len(angledict),4],dtype='float')

  for a in range(len(angledict)):
    temp=angledict[a+1]
    astr1=temp[0:3]
    lmptype=int(a+1)
    for c in range(len(ANGLE_ARRAY)):
      astr2=[atom_type[ANGLE_ARRAY[c,0]-1][0],atom_type[ANGLE_ARRAY[c,1]-1][0],atom_type[ANGLE_ARRAY[c,2]-1][0]]
      astr2R=astr2[::-1]
      if cmp(astr1,astr2)==0:
        angletypes.append(int(lmptype))

        x1_a=coord_regblah[ANGLE_ARRAY[c,0]-1,0]
        x2_a=coord_regblah[ANGLE_ARRAY[c,1]-1,0]
        x3_a=coord_regblah[ANGLE_ARRAY[c,2]-1,0]
        y1_a=coord_regblah[ANGLE_ARRAY[c,0]-1,1]
        y2_a=coord_regblah[ANGLE_ARRAY[c,1]-1,1]
        y3_a=coord_regblah[ANGLE_ARRAY[c,2]-1,1]
        z1_a=coord_regblah[ANGLE_ARRAY[c,0]-1,2]
        z2_a=coord_regblah[ANGLE_ARRAY[c,1]-1,2]
        z3_a=coord_regblah[ANGLE_ARRAY[c,2]-1,2]

        dist1=np.sqrt((x1_a-x2_a)**2+(y1_a-y2_a)**2+(z1_a-z2_a)**2)
        dist2=np.sqrt((x1_a-x3_a)**2+(y1_a-y3_a)**2+(z1_a-z3_a)**2)
        dist3=np.sqrt((x3_a-x2_a)**2+(y3_a-y2_a)**2+(z3_a-z2_a)**2)

        if dist1 <= 10 and dist2 <= 10 and dist3 <= 10:
          vector_one=(x1_a - x2_a,y1_a - y2_a,z1_a - z2_a)
          vector_two=(x3_a - x2_a,y3_a - y2_a,z3_a - z2_a)
          angle=np.degrees(angle_between(vector_one,vector_two))
          fileANGLES.write(str(lmptype)+' '+str(angle)+'\n')

          list_angle_atom_types=angledict[lmptype]
          angle_infoA[lmptype-1,0]=list_angle_atom_types[0]
          angle_infoA[lmptype-1,1]=list_angle_atom_types[1]
          angle_infoA[lmptype-1,2]=list_angle_atom_types[2]
          angle_infoB[lmptype-1,0]=lmptype
          angle_infoB[lmptype-1,1]=angle_infoB[lmptype-1,1]+angle
          angle_infoB[lmptype-1,2]=angle_infoB[lmptype-1,2]+1
          angle_infoB[lmptype-1,3]=float(angle_infoB[lmptype-1,1])/float(angle_infoB[lmptype-1,2])

        REALANGLEARRAY=np.vstack((REALANGLEARRAY,[ANGLE_ARRAY[c,0],ANGLE_ARRAY[c,1],ANGLE_ARRAY[c,2]]))
      elif cmp(astr1,astr2R)==0:
        angletypes.append(int(lmptype))

        x1_a=coord_regblah[ANGLE_ARRAY[c,0]-1,0]
        x2_a=coord_regblah[ANGLE_ARRAY[c,1]-1,0]
        x3_a=coord_regblah[ANGLE_ARRAY[c,2]-1,0]
        y1_a=coord_regblah[ANGLE_ARRAY[c,0]-1,1]
        y2_a=coord_regblah[ANGLE_ARRAY[c,1]-1,1]
        y3_a=coord_regblah[ANGLE_ARRAY[c,2]-1,1]
        z1_a=coord_regblah[ANGLE_ARRAY[c,0]-1,2]
        z2_a=coord_regblah[ANGLE_ARRAY[c,1]-1,2]
        z3_a=coord_regblah[ANGLE_ARRAY[c,2]-1,2]

        dist1=np.sqrt((x1_a-x2_a)**2+(y1_a-y2_a)**2+(z1_a-z2_a)**2)
        dist2=np.sqrt((x1_a-x3_a)**2+(y1_a-y3_a)**2+(z1_a-z3_a)**2)
        dist3=np.sqrt((x3_a-x2_a)**2+(y3_a-y2_a)**2+(z3_a-z2_a)**2)

        if dist1 <= 10 and dist2 <= 10 and dist3 <= 10:
          vector_one=(x1_a - x2_a,y1_a - y2_a,z1_a - z2_a)
          vector_two=(x3_a - x2_a,y3_a - y2_a,z3_a - z2_a)
          angle=np.degrees(angle_between(vector_one,vector_two))
          fileANGLES.write(str(lmptype)+' '+str(angle)+'\n')

          list_angle_atom_types=angledict[lmptype]
          angle_infoA[lmptype-1,0]=list_angle_atom_types[0]
          angle_infoA[lmptype-1,1]=list_angle_atom_types[1]
          angle_infoA[lmptype-1,2]=list_angle_atom_types[2]
          angle_infoB[lmptype-1,0]=lmptype
          angle_infoB[lmptype-1,1]=angle_infoB[lmptype-1,1]+angle
          angle_infoB[lmptype-1,2]=angle_infoB[lmptype-1,2]+1
          angle_infoB[lmptype-1,3]=float(angle_infoB[lmptype-1,1])/float(angle_infoB[lmptype-1,2])

        REALANGLEARRAY=np.vstack((REALANGLEARRAY,[ANGLE_ARRAY[c,0],ANGLE_ARRAY[c,1],ANGLE_ARRAY[c,2]]))
  print('REAL NUMBER OF ANGLES: '+str(len(angletypes)))
  print('SIZE OF CORRECT ANGLE ARRAY: '+str(len(REALANGLEARRAY)))

  fileANGLES.close()

  fileANGLES=open(output+'/SUMMARY_ANGLES_TYPED_LENGTHS.data','w')
  for b in range(len(angledict)):
    fileANGLES.write(str(angle_infoA[b,0])+' '+str(angle_infoA[b,1])+' '+str(angle_infoA[b,2])+' '+str(int(angle_infoB[b,0]))+' '+str(angle_infoB[b,1])+' '+str(int(angle_infoB[b,2]))+' '+str(angle_infoB[b,3])+'\n')      
  fileANGLES.close()

  return(angletypes,REALANGLEARRAY)

#DIHEDRAL TYPING
def dihedraltyping(dihedraldict,DIHEDRAL_ARRAY,atom_type):
  REALDIHEDRALARRAY=np.ones([0,4],dtype='int')

  DIHED_FLAGS=np.zeros([len(DIHEDRAL_ARRAY),1],dtype='int')

  dihedraltypes=[]
  for d in range(len(dihedraldict)):
    temp=dihedraldict[d+1]
    astr1=temp[0:4]
    lmptype=d+1
    for c in range(len(DIHEDRAL_ARRAY)):
      astr2=[atom_type[DIHEDRAL_ARRAY[c,0]-1][0],atom_type[DIHEDRAL_ARRAY[c,1]-1][0],
             atom_type[DIHEDRAL_ARRAY[c,2]-1][0],atom_type[DIHEDRAL_ARRAY[c,3]-1][0]]
      astr2R=astr2[::-1]
      astr3=[atom_type[DIHEDRAL_ARRAY[c,0]-1][0],atom_type[DIHEDRAL_ARRAY[c,2]-1][0],
             atom_type[DIHEDRAL_ARRAY[c,1]-1][0],atom_type[DIHEDRAL_ARRAY[c,3]-1][0]]
      astr3R=astr3[::-1]
      #4 TOTAL possible ways to write the same angle
      if cmp(astr1,astr2)==0:
        dihedraltypes.append(int(lmptype))
        REALDIHEDRALARRAY=np.vstack((REALDIHEDRALARRAY,[DIHEDRAL_ARRAY[c,0],DIHEDRAL_ARRAY[c,1],DIHEDRAL_ARRAY[c,2],DIHEDRAL_ARRAY[c,3]]))
        DIHED_FLAGS[c,0]=1
      elif cmp(astr1,astr2R)==0:
        dihedraltypes.append(int(lmptype))
        REALDIHEDRALARRAY=np.vstack((REALDIHEDRALARRAY,[DIHEDRAL_ARRAY[c,0],DIHEDRAL_ARRAY[c,1],DIHEDRAL_ARRAY[c,2],DIHEDRAL_ARRAY[c,3]]))
        DIHED_FLAGS[c,0]=1
      elif cmp(astr1,astr3)==0:
        dihedraltypes.append(int(lmptype))
        REALDIHEDRALARRAY=np.vstack((REALDIHEDRALARRAY,[DIHEDRAL_ARRAY[c,0],DIHEDRAL_ARRAY[c,1],DIHEDRAL_ARRAY[c,2],DIHEDRAL_ARRAY[c,3]]))
        DIHED_FLAGS[c,0]=1
      elif cmp(astr1,astr3R)==0:
        dihedraltypes.append(int(lmptype))
        REALDIHEDRALARRAY=np.vstack((REALDIHEDRALARRAY,[DIHEDRAL_ARRAY[c,0],DIHEDRAL_ARRAY[c,1],DIHEDRAL_ARRAY[c,2],DIHEDRAL_ARRAY[c,3]]))
        DIHED_FLAGS[c,0]=1

  print('REAL NUMBER OF PROPER DIHEDRALS: '+str(len(dihedraltypes)))
  print('SIZE OF CORRECT DIHEDRAL ARRAY: '+str(len(REALDIHEDRALARRAY)))
  print('SIZE OF CORRECT DIHEDRAL ARRAY2: '+str(np.sum(DIHED_FLAGS)))

  for c in range(len(DIHEDRAL_ARRAY)):
    if DIHED_FLAGS[c,0]==0:
      astr2=[atom_type[DIHEDRAL_ARRAY[c,0]-1][0],atom_type[DIHEDRAL_ARRAY[c,1]-1][0],
             atom_type[DIHEDRAL_ARRAY[c,2]-1][0],atom_type[DIHEDRAL_ARRAY[c,3]-1][0]]
      print(astr2)

  return(dihedraltypes,REALDIHEDRALARRAY)

#IMPROPER TYPING
def impropertyping(improperdict,atom_type,conn_matrix_char_TYPED_HYBRID,conn_matrix_nums,centeratomINDEX):
  REALIMPROPERARRAY=np.ones([0,4],dtype='int')
  impropertypes=[]

  for improper in range(len(improperdict)):
    temp=improperdict[improper+1]
    astr1=temp[0:4]
    if improper_proper_flag=='YES-improper-as-proper':
      lmptype=improper+1+ntorsiontypes
    centerimproperatom=(temp[centeratomINDEX-1])

    for atom in range(len(atom_type)):
      if atom_type[atom][0] == centerimproperatom: #check to see if the atom is one of the central atoms 
        #Find all its nearest neighbors from the conn_matrix_char and conn_matrix_nums
        #Place in the improper array in the correct positions
        if centeratomINDEX == 3:
          astr2=[conn_matrix_char_TYPED_HYBRID[atom][0],conn_matrix_char_TYPED_HYBRID[atom][1],atom_type[atom][0],conn_matrix_char_TYPED_HYBRID[atom][2]]
          astr3=[conn_matrix_char_TYPED_HYBRID[atom][1],conn_matrix_char_TYPED_HYBRID[atom][0],atom_type[atom][0],conn_matrix_char_TYPED_HYBRID[atom][2]]

          astr4=[conn_matrix_char_TYPED_HYBRID[atom][0],conn_matrix_char_TYPED_HYBRID[atom][2],atom_type[atom][0],conn_matrix_char_TYPED_HYBRID[atom][1]]
          astr5=[conn_matrix_char_TYPED_HYBRID[atom][2],conn_matrix_char_TYPED_HYBRID[atom][0],atom_type[atom][0],conn_matrix_char_TYPED_HYBRID[atom][1]]

          astr6=[conn_matrix_char_TYPED_HYBRID[atom][1],conn_matrix_char_TYPED_HYBRID[atom][2],atom_type[atom][0],conn_matrix_char_TYPED_HYBRID[atom][0]]
          astr7=[conn_matrix_char_TYPED_HYBRID[atom][2],conn_matrix_char_TYPED_HYBRID[atom][1],atom_type[atom][0],conn_matrix_char_TYPED_HYBRID[atom][0]]

          if list(astr1)==list(astr2):
            REALIMPROPERARRAY=np.vstack((REALIMPROPERARRAY,[conn_matrix_nums[atom][0],conn_matrix_nums[atom][1],atom+1,conn_matrix_nums[atom][2]]))
            impropertypes.append(int(lmptype))
          elif list(astr1)==list(astr3):
            REALIMPROPERARRAY=np.vstack((REALIMPROPERARRAY,[conn_matrix_nums[atom][1],conn_matrix_nums[atom][0],atom+1,conn_matrix_nums[atom][2]]))
            impropertypes.append(int(lmptype))
          elif list(astr1)==list(astr4):
            REALIMPROPERARRAY=np.vstack((REALIMPROPERARRAY,[conn_matrix_nums[atom][0],conn_matrix_nums[atom][2],atom+1,conn_matrix_nums[atom][1]]))
            impropertypes.append(int(lmptype))
          elif list(astr1)==list(astr5):
            REALIMPROPERARRAY=np.vstack((REALIMPROPERARRAY,[conn_matrix_nums[atom][2],conn_matrix_nums[atom][0],atom+1,conn_matrix_nums[atom][1]]))
            impropertypes.append(int(lmptype))
          elif list(astr1)==list(astr6):
            REALIMPROPERARRAY=np.vstack((REALIMPROPERARRAY,[conn_matrix_nums[atom][1],conn_matrix_nums[atom][2],atom+1,conn_matrix_nums[atom][0]]))
            impropertypes.append(int(lmptype))
          elif list(astr1)==list(astr7):
            REALIMPROPERARRAY=np.vstack((REALIMPROPERARRAY,[conn_matrix_nums[atom][2],conn_matrix_nums[atom][1],atom+1,conn_matrix_nums[atom][0]]))
            impropertypes.append(int(lmptype))

  print('REAL NUMBER OF IMPROPER DIHEDRALS: '+str(len(impropertypes)))
  print('SIZE OF CORRECT IMPROPER ARRAY: '+str(len(REALIMPROPERARRAY)))
  return(impropertypes,REALIMPROPERARRAY)

def lammps_write(myid,name,natoms,atom_type,LAMMPS_ATOMINDEXESblah,LAMMPS_CHARGESblah,coord_regblah,BOND_ARRAY,ANGLE_ARRAY,DIHEDRAL_ARRAY,checkflag,conn_matrix_char_TYPED_HYBRID,conn_matrix_nums):
    file1=open(output+'/'+str(name)+'.data','w')

    bondtypes,REALBONDARRAY=bondtyping(bonddict,BOND_ARRAY,atom_type,coord_regblah)
    angletypes,REALANGLEARRAY=angletyping(angledict,ANGLE_ARRAY,atom_type,coord_regblah)
    dihedraltypes,REALDIHEDRALARRAY=dihedraltyping(dihedraldict,DIHEDRAL_ARRAY,atom_type)
    impropertypes,REALIMPROPERARRAY=impropertyping(improperdict,atom_type,conn_matrix_char_TYPED_HYBRID,conn_matrix_nums,centeratomINDEX)
    

    file1.write('LAMMPS DATA FILE FOR '+str(name)+'\n\n')

    file1.write(str(natoms)+' ' +'atoms\n')
    file1.write(str(len(REALBONDARRAY))+' ' +'bonds\n')
    file1.write(str(len(REALANGLEARRAY))+' ' +'angles\n')
    if improper_proper_flag=='YES-improper-as-proper':
      file1.write(str(len(REALDIHEDRALARRAY)+len(REALIMPROPERARRAY))+' ' +'dihedrals\n\n')
    else:
      file1.write(str(len(REALDIHEDRALARRAY))+' ' +'dihedrals\n\n')

    file1.write(str(natypes)+' '+'atom types\n')
    file1.write(str(nbondtypes)+' ' +'bond types\n')
    file1.write(str(nangletypes)+' ' +'angle types\n')
    if improper_proper_flag=='YES-improper-as-proper':
      file1.write(str(ntorsiontypes+nimpropertypes)+' ' +'dihedral types\n\n')
    else:
      file1.write(str(ntorsiontypes)+' ' +'dihedral types\n\n')

    file1.write(str(0.0)+' '+str(lx)+' '+'xlo xhi\n')
    file1.write(str(0.0)+' '+str(ly)+' '+'ylo yhi\n')
    file1.write(str(0.0)+' '+str(lz)+' '+'zlo zhi\n')
    if triclinic_flag=='YES':
      file1.write(str(xytilt)+' '+str(xztilt)+' '+str(yztilt)+' xy xz yz\n\n')
    if triclinic_flag=='NO':
      file1.write('\n')

    file1.write('Masses\n\n')
    for atom in range(natypes):
        file1.write(str(atom+1)+' '+str(masses[atom])+'\n')

    file1.write('\n')
    file1.write('Bond Coeffs\n\n')
    for bond in range(nbondtypes):
        file1.write(str(bond+1)+' '+str(float(bond_param[bond,0])/float(bond_KCORRECTION))+' '+str(bond_param[bond,1])+'\n') #------------------------> INPUT FILE

    file1.write('\n')
    file1.write('Angle Coeffs\n\n')
    for angle in range(nangletypes):
        file1.write(str(angle+1)+' '+str(float(angle_param[angle,0])/float(angle_KCORRECTION))+' '+str(angle_param[angle,1])+'\n') #------------------------> INPUT FILE

    file1.write('\n')
    file1.write('Dihedral Coeffs\n\n')
    for dihedral in range(ntorsiontypes):
        file1.write(str(dihedral+1)+' '+str(torsion_param[dihedral,0])+' '+str(int(torsion_param[dihedral,2]))
                                  +' '+str(int(torsion_param[dihedral,1]))+' '+str(0)+'\n')

    if improper_proper_flag=='YES-improper-as-proper':
      for improper in range(nimpropertypes):
        file1.write(str(dihedral+improper+2)+' '+str(improper_param[improper,0])+' '+str(int(improper_param[improper,2]))
                                   +' '+str(int(improper_param[improper,1]))+' '+str(0)+'\n')


    #ATOMS
    if checkflag=='YES-checktyping':
      file1.write('\n')
      file1.write('Atoms\n\n')
      for atom in range(natoms):
          file1.write(str(atom+1)+' '+str(1)+' '+str(atom_type[atom][0])+' '+
                     str(LAMMPS_ATOMINDEXESblah[atom]).replace("[",'').replace("]",'')+' '+
                     str(LAMMPS_CHARGESblah[atom]).replace("[",'').replace("]",'')+' '+
                     str(coord_regblah[atom,0])+' '+str(coord_regblah[atom,1])+' '+str(coord_regblah[atom,2])+' '+
                     str(0)+' '+str(0)+' '+str(0)+' '+'\n')
    else:
      file1.write('\n')
      file1.write('Atoms\n\n')
      for atom in range(natoms):
          file1.write(str(atom+1)+' '+str(1)+' '+
                     str(LAMMPS_ATOMINDEXESblah[atom]).replace("[",'').replace("]",'')+' '+
                     str(LAMMPS_CHARGESblah[atom]).replace("[",'').replace("]",'')+' '+
                     str(coord_regblah[atom,0])+' '+str(coord_regblah[atom,1])+' '+str(coord_regblah[atom,2])+' '+
                     str(0)+' '+str(0)+' '+str(0)+' '+'\n')

    #VELOCITIES
    file1.write('\n')
    file1.write('Velocities\n\n')
    for veloc in range(natoms):
        file1.write(str(veloc+1)+' '+str(0)+' '+str(0)+' '+str(0)+'\n')

    #BONDS
    if checkflag=='YES-checktyping':
      file1.write('\n')
      file1.write('Bonds\n\n')
      for bond in range(len(bondtypes)):
          file1.write(str(bond+1)+' '+str(bondtypes[bond])+' '+str(REALBONDARRAY[bond][0])+' '+
                                                               str(REALBONDARRAY[bond][1])+' '+
                                                               str(atom_type[REALBONDARRAY[bond][0]-1][0])+' '+
                                                               str(atom_type[REALBONDARRAY[bond][1]-1][0])+'\n')
    else:
      file1.write('\n')
      file1.write('Bonds\n\n')
      for bond in range(len(bondtypes)):
          file1.write(str(bond+1)+' '+str(bondtypes[bond])+' '+str(REALBONDARRAY[bond][0])+' '+
                                                               str(REALBONDARRAY[bond][1])+'\n')

    #ANGLES
    if checkflag=='YES-checktyping':
      file1.write('\n')
      file1.write('Angles\n\n')
      for angle in range(len(angletypes)):
          file1.write(str(angle+1)+' '+str(angletypes[angle])+' '+str(REALANGLEARRAY[angle][0])+' '+
                                                                  str(REALANGLEARRAY[angle][1])+' '+
                                                                  str(REALANGLEARRAY[angle][2])+' '+
                                                                  str(atom_type[REALANGLEARRAY[angle][0]-1][0])+' '+
                                                                  str(atom_type[REALANGLEARRAY[angle][1]-1][0])+' '+
                                                                  str(atom_type[REALANGLEARRAY[angle][2]-1][0])+'\n')
    else:
      file1.write('\n')
      file1.write('Angles\n\n')
      for angle in range(len(angletypes)):
          file1.write(str(angle+1)+' '+str(angletypes[angle])+' '+str(REALANGLEARRAY[angle][0])+' '+
                                                                  str(REALANGLEARRAY[angle][1])+' '+
                                                                  str(REALANGLEARRAY[angle][2])+'\n')

    #PROPERS
    if checkflag=='YES-checktyping':
      file1.write('\n')
      file1.write('Dihedrals\n\n')
      for proper in range(len(dihedraltypes)):
          file1.write(str(proper+1)+' '+str(dihedraltypes[proper])+' '+str(REALDIHEDRALARRAY[proper][0])+' '+
                                                                       str(REALDIHEDRALARRAY[proper][1])+' '+
                                                                       str(REALDIHEDRALARRAY[proper][2])+' '+
                                                                       str(REALDIHEDRALARRAY[proper][3])+' '+
                                                                       str(atom_type[REALDIHEDRALARRAY[proper][0]-1][0])+' '+
                                                                       str(atom_type[REALDIHEDRALARRAY[proper][1]-1][0])+' '+
                                                                       str(atom_type[REALDIHEDRALARRAY[proper][2]-1][0])+' '+
                                                                       str(atom_type[REALDIHEDRALARRAY[proper][3]-1][0])+'\n')    
    else:
      file1.write('\n')
      file1.write('Dihedrals\n\n')
      for proper in range(len(dihedraltypes)):
          file1.write(str(proper+1)+' '+str(dihedraltypes[proper])+' '+str(REALDIHEDRALARRAY[proper][0])+' '+
                                                                       str(REALDIHEDRALARRAY[proper][1])+' '+
                                                                       str(REALDIHEDRALARRAY[proper][2])+' '+
                                                                       str(REALDIHEDRALARRAY[proper][3])+'\n')

    #IMPROPERS
    if improper_proper_flag=='YES-improper-as-proper':
      if checkflag=='YES-checktyping':
        for improper in range(len(impropertypes)):
            file1.write(str(proper+improper+2)+' '+str(impropertypes[improper])+' '+str(REALIMPROPERARRAY[improper][0])+' '+
                                                                                    str(REALIMPROPERARRAY[improper][1])+' '+
                                                                                    str(REALIMPROPERARRAY[improper][2])+' '+
                                                                                    str(REALIMPROPERARRAY[improper][3])+' '+
                                                                                    str(atom_type[REALIMPROPERARRAY[improper][0]-1][0])+' '+
                                                                                    str(atom_type[REALIMPROPERARRAY[improper][1]-1][0])+' '+
                                                                                    str(atom_type[REALIMPROPERARRAY[improper][2]-1][0])+' '+
                                                                                    str(atom_type[REALIMPROPERARRAY[improper][3]-1][0])+'\n')    
      else:
        for improper in range(len(impropertypes)):
            file1.write(str(proper+improper+2)+' '+str(impropertypes[improper])+' '+str(REALIMPROPERARRAY[improper][0])+' '+
                                                                                    str(REALIMPROPERARRAY[improper][1])+' '+
                                                                                    str(REALIMPROPERARRAY[improper][2])+' '+
                                                                                    str(REALIMPROPERARRAY[improper][3])+'\n')


    file1.close()
    return()

############################################################################################################
def cif_write(myid,name,natoms,ATOM_TYPE_REAL,coord_sorted_frac,aLEN,bLEN,cLEN,alphaD,betaD,gammaD,x_expand,y_expand,z_expand):
    #Writing the CIF file
    file=open(output+'/'+str(name)+'.cif','w')
    file.write('#======================================================================\n\n')
    file.write('# CRYSTAL DATA\n\n')
    file.write('#----------------------------------------------------------------------\n\n')
    file.write('data_VESTA_phase_1\n\n\n')
    file.write('_pd_phase_name                         '+str(name)+'\n')
    file.write('_cell_length_a                         '+str(aLEN*x_expand)+'\n')
    file.write('_cell_length_b                         '+str(bLEN*y_expand)+'\n')
    file.write('_cell_length_c                         '+str(cLEN*z_expand)+'\n')
    file.write('_cell_angle_alpha                      '+str(alphaD)+'\n')
    file.write('_cell_angle_beta                       '+str(betaD)+'\n')
    file.write('_cell_angle_gamma                      '+str(gammaD)+'\n')
    file.write('_symmetry_space_group_name_H-M         "P 1"\n')
    file.write('_symmetry_Int_Tables_number            1\n\n')
    file.write('loop_\n')
    file.write('   "x, y, z"\n\n')
    file.write('loop_\n')
    file.write('   _atom_site_label\n')
    file.write('   _atom_site_occupancy\n')
    file.write('   _atom_site_fract_x\n')
    file.write('   _atom_site_fract_y\n')
    file.write('   _atom_site_fract_z\n')
    file.write('   _atom_site_B_iso_or_equiv\n')
    file.write('   _atom_site_adp_type\n')
    file.write('   _atom_site_type_symbol\n')
    for i in range(natoms):
        file.write(str(ATOM_TYPE_REAL[i]).replace("['",'').replace("']",'')+' 1.0 '+str(coord_sorted_frac[i,0])+' '+str(coord_sorted_frac[i,1])+' '+str(coord_sorted_frac[i,2])+
                   ' Biso 1.0 '+str(ATOM_TYPE_REAL[i]).replace("['",'').replace("']",'')+'\n')
    file.write('\r')
    file.close()
    print('CIF file has been written.')
    return(1)

def cif_writeALTERNATE(myid,name,natoms,ATOM_TYPE_REAL,coord_sorted_frac,aLEN,bLEN,cLEN,alphaD,betaD,gammaD,x_expand,y_expand,z_expand):
    #Writing the CIF file
    file=open(output+'/'+str(name)+'.cif','w')
    file.write('#======================================================================\n\n')
    file.write('# CRYSTAL DATA\n\n')
    file.write('#----------------------------------------------------------------------\n\n')
    file.write('data_VESTA_phase_1\n\n\n')
    file.write('_pd_phase_name                         '+str(name)+'\n')
    file.write('_cell_length_a                         '+str(aLEN*x_expand)+'\n')
    file.write('_cell_length_b                         '+str(bLEN*y_expand)+'\n')
    file.write('_cell_length_c                         '+str(cLEN*z_expand)+'\n')
    file.write('_cell_angle_alpha                      '+str(alphaD)+'\n')
    file.write('_cell_angle_beta                       '+str(betaD)+'\n')
    file.write('_cell_angle_gamma                      '+str(gammaD)+'\n')
    file.write('_symmetry_space_group_name_H-M         "P 1"\n')
    file.write('_symmetry_Int_Tables_number            1\n\n')
    file.write('loop_\n')
    file.write('   "x, y, z"\n\n')
    file.write('loop_\n')
    file.write('   _atom_site_label\n')
    file.write('   _atom_site_occupancy\n')
    file.write('   _atom_site_fract_x\n')
    file.write('   _atom_site_fract_y\n')
    file.write('   _atom_site_fract_z\n')
    file.write('   _atom_site_B_iso_or_equiv\n')
    file.write('   _atom_site_adp_type\n')
    file.write('   _atom_site_type_symbol\n')
    for i in range(natoms):
        file.write(str(ATOM_TYPE_REAL[i]).replace("['",'').replace("']",'')+' 1.0 '+str(coord_sorted_frac[0,i])+' '+str(coord_sorted_frac[1,i])+' '+str(coord_sorted_frac[2,i])+
                   ' Biso 1.0 '+str(ATOM_TYPE_REAL[i]).replace("['",'').replace("']",'')+'\n')
    file.write('\r')
    file.close()
    print('CIF file has been written.')
    return(1)
############################################################################################################
def measure_SOP(linker_composition,linker_binaryNNlist,n_linkerNNs):
    SOP_store=[]
    for clust in range(len(linker_binaryNNlist)):
        if linker_binaryNNlist[clust,0]==0: #WHY WOULD I STAND ON ZIF-90 LINKERS????????????????
            probAB=float(np.sum(linker_binaryNNlist[clust,1:n_linkerNNs+1]))/float(n_linkerNNs) #THIS IS THE PROBABILITY THAT I SEE A ZIF-90 (PD)
            if linker_composition != float(0):
                SOP=1-float(probAB)/float(linker_composition)
            else:
                SOP=0
            SOP_store.append(SOP)
    SOP_store=np.array(SOP_store,dtype='float')
    avg_SOP=np.average(SOP_store)
    return(avg_SOP)

#def switch_linkers(rand_linker_1,rand_linker_2,NCL,linker_NNlist,linkersswitching,n_linkerNNs): #THIS NEEDS TO BE PLACED INTO FORTRAN
#    #Replace the occurances of the linker 2 with linker 1 (it now becomes a 90)
#    LSlist=[rand_linker_1 if x==rand_linker_2 else x for x in linkersswitching]

#    #Do the same for each row of the linker_binaryNNlist
#    LBNNlist=np.zeros((NCL,n_linkerNNs+1),dtype='int')
#    for clust in range(NCL):
#         clustID=clust+1
#         if clustID in LSlist:
#             LBNNlist[clust,0]=1
#         for clustnum in range(NCL):
#             clustIDB=clustnum+1
#             for neighbor in range(n_linkerNNs):
#                 if clustIDB==linker_NNlist[clust,neighbor] and (clustIDB in LSlist):
#                     LBNNlist[clust,neighbor+1]=1
#    return(LBNNlist,LSlist)

def switch_linkers(rand_linker_1,rand_linker_2,NCL,linker_NNlist,linkersswitching,n_linkerNNs,LBNNlistB): #THIS NEEDS TO BE PLACED INTO FORTRAN
    #Replace the occurances of the linker 2 with linker 1 (it now becomes a 90)
    #print(rand_linker_1,rand_linker_2)
    #print(linkersswitching)
    LSlist=[rand_linker_1 if x==rand_linker_2 else x for x in linkersswitching]
    #print(LSlist)

    LBNNlist=copy.deepcopy(LBNNlistB)
    NNs_linkerswitch_1=linker_NNlist[rand_linker_1-1,:] #these get turned to 1's
    NNs_linkerswitch_2=linker_NNlist[rand_linker_2-1,:] #these get turned to 0's
    
    #print(rand_linker_1,rand_linker_2)
    #print(NNs_linkerswitch_1)
    #print(NNs_linkerswitch_2)


    #turning to 1's
    LBNNlist[rand_linker_1-1,0]=1
    for i in range(len(NNs_linkerswitch_1)): #this is always 6
      NNindex=int(NNs_linkerswitch_1[i])-1
      for neighbor in range(n_linkerNNs): #this is always 6
        if rand_linker_1==linker_NNlist[NNindex,neighbor]:
          LBNNlist[NNindex,neighbor+1]=1

    #turning to 0's
    LBNNlist[rand_linker_2-1,0]=0
    for i in range(len(NNs_linkerswitch_2)): #this is always 6
      NNindex=int(NNs_linkerswitch_2[i])-1
      for neighbor in range(n_linkerNNs): #this is always 6
        if rand_linker_2==linker_NNlist[NNindex,neighbor]:
          LBNNlist[NNindex,neighbor+1]=0

    #Do the same for each row of the linker_binaryNNlist
    #LBNNlist=np.zeros((NCL,n_linkerNNs+1),dtype='int')
    #for clust in range(NCL):
    #     clustID=clust+1
    #     if clustID in LSlist:
    #         LBNNlist[clust,0]=1
    #     for clustnum in range(NCL):
    #         clustIDB=clustnum+1
    #         for neighbor in range(n_linkerNNs):
    #             if clustIDB==linker_NNlist[clust,neighbor] and (clustIDB in LSlist):
    #                 LBNNlist[clust,neighbor+1]=1
    return(LBNNlist,LSlist)
############################################################################################################


def main(myid,main_argument):
  ###Set the main paths 
     global running_dir
     running_dir=str(os.getcwd()).replace('/src','')
     global input
     input=running_dir+'/input'
     global output
     output=running_dir+'/output'
     global src
     src=running_dir+'/src' 
     global databases
     databases=running_dir+'/databases'


  ###Check to see if the appropriate folders are available
     #if input:
          #os.mkdir()
     #if output:
          #os.mkdir()
     #if databases:
          #os.mkdir()

  ###Read input files
     sigma,set_SOP=read_inputs()
     #print(globals())
     #print(locals())
     file_name=input+'/'+str(file_argument)
     name=file_argument.replace('.xyz','')   

  ###Set the random seed
     np.random.seed(GLOBAL_SEED)#--------------------------------------------------------------------------->DO IN INPUT SCRIPT  

  ###Read in the databases (color, VDW radii, molecular weight, periodic table location, etc.)
     oneprint(myid,'Database read.')
     atype,radii=database_read(databases)

  ###Read xyz file
     oneprint(myid,'XYZ read.')
     print(file_name)
     atom_type,xpos,ypos,zpos=xyz_read(file_name) #do file format checking to catch errors
     natoms=len(atom_type)

  ###Expand the unit cell
     xposnpc=np.array(xpos,dtype='f')
     yposnpc=np.array(ypos,dtype='f')
     zposnpc=np.array(zpos,dtype='f')
     coord_reg=np.vstack([xposnpc,yposnpc,zposnpc])
     coord_reg=coord_reg.T
     print(coord_reg.shape)

     C_mat,C_mat_inv=transformation_matrix(aLEN,bLEN,cLEN,alphaD,betaD,gammaD,1,1,1)
     C_mat_LARGE,C_mat_inv_LARGE=transformation_matrix(aLEN,bLEN,cLEN,alphaD,betaD,gammaD,x_expand,y_expand,z_expand)
     coord_frac=np.dot(C_mat_inv,coord_reg.T)
     coord_frac=np.array(coord_frac.T)
     print(coord_frac.shape)

     atom_type,coord_frac=expand_UC(myid,natoms,atom_type,coord_frac,x_expand,y_expand,z_expand)
     coord_frac=np.array(coord_frac.T)
     print(coord_frac.shape)
     coord_reg=np.dot(C_mat,coord_frac.T)
     print(coord_reg.shape)
     coord_reg=coord_reg.T
     print(coord_reg.shape)
     xpos=coord_reg[:,0]
     ypos=coord_reg[:,1]
     zpos=coord_reg[:,2]

     atom_type=np.array(atom_type,dtype='S7')
     natoms=len(atom_type)
     atom_typeINT=[]
     for atom in range(natoms):
         if atom_type[atom]=='H' or atom_type[atom]=='H2' or atom_type[atom]=='H3':
             atom_typeINT.append(1)
         if atom_type[atom]=='C' or atom_type[atom]=='C1' or atom_type[atom]=='C2' or atom_type[atom]=='C3':
             atom_typeINT.append(2)
         if atom_type[atom]=='N':
             atom_typeINT.append(3)
         if atom_type[atom]=='O':
             atom_typeINT.append(4)
         if atom_type[atom]=='Zn':
             atom_typeINT.append(5)
     atom_typeINT=np.array(atom_typeINT,dtype='int')

     file=open(output+'/'+str(name)+'TEST_EXPAND'+str(x_expand)+str(y_expand)+str(z_expand)+'.xyz','w')
     file.write(str(len(atom_type))+'\n')
     file.write(str(name)+' '+str(aLEN)+' '+str(bLEN)+' '+str(cLEN)+' '+str(alphaD)+' '+str(betaD)+' '+str(gammaD)+'\n')
     for i in range(len(atom_type)):
         file.write(str(atom_type[i]).replace("['",'').replace("']",'')+' '+str(float(xpos[i]))+' '+str(float(ypos[i]))+' '+str(float(zpos[i]))+'\n')
     file.write('\r')
     file.close()

     file=open(output+'/TEST_FRAC'+str(x_expand)+str(y_expand)+str(z_expand)+'.xyz','w')
     file.write(str(natoms)+'\n\n')
     for x in range(natoms):
         #file.write(str(atom_type[x])+' '+str(coord_frac[x,0]*aLEN)+' '+str(coord_frac[x,1]*bLEN)+' '+str(coord_frac[x,2]*cLEN)+'\n')
         file.write(str(atom_type[x])+' '+str(coord_frac[x,0]*1)+' '+str(coord_frac[x,1]*1)+' '+str(coord_frac[x,2]*1)+'\n')
     file.close()

  ###Remove duplicate atoms 

  ###Generate the connectivity matrix and bonds
     oneprint(myid,'Connectivity Map generating.')
     #this will need the expansion in xyz directions when I expand unit cell
     conn_matrix1=np.zeros([natoms,13])
     conn_matrix_char=np.zeros([natoms,6],dtype='a3')
     conn_matrix_nums1=np.zeros([natoms,6],dtype='int')
     bond_fix = bond_corr
     radii_array=[]
     for atom in range(natoms):
         atom1=atom_type[atom]
         rad1=cov_radius_assignment(atype,radii,atom1)
         radii_array.append(rad1)
     radii_array1=np.array(radii_array,dtype='f')
     print (COMBINEDv2.connectivity_mapping.__doc__)
     conn_matrix1,conn_matrix_nums1=COMBINEDv2.connectivity_mapping(atom_typeINT,radii_array1,coord_frac,C_mat,x_expand,y_expand,z_expand,conn_matrix1,conn_matrix_nums1,bond_fix)
     #conn_matrix,conn_matrix_nums,conn_matrix_char=connectivity_mapping(natoms,radii,atom_type,atype,coord_frac,C_mat,x_expand,y_expand,z_expand) 
     #print(conn_matrix_char)
     #print(conn_matrix_nums1)
     #print(conn_matrix1)
     conn_matrix=copy.deepcopy(conn_matrix1)
     conn_matrix_nums=copy.deepcopy(conn_matrix_nums1)

     file=open(output+'/connect_check2_nums.xyz','w')
     for x in range(natoms):
         file.write(str(conn_matrix_nums[x,0])+' '+str(conn_matrix_nums[x,1])+' '+str(conn_matrix_nums[x,2])+' '+str(conn_matrix_nums[x,3])+' '+str(conn_matrix_nums[x,4])+' '+str(conn_matrix_nums[x,5])+'\n')
     file.close()

     for atom in range(natoms):
         for i in range(6):
             if conn_matrix_nums[atom,i] != 0 and conn_matrix_nums[atom,i] <= natoms:
                 conn_matrix_char[atom,i]=atom_type[int(conn_matrix_nums[atom,i])-1]
     #print(conn_matrix_char)

     oneprint(myid,'Connectivity Map generated!')

     file=open(output+'/connect_check1_char.xyz','w')
     for x in range(natoms):
         file.write(str(conn_matrix_char[x,0])+' '+str(conn_matrix_char[x,1])+' '+str(conn_matrix_char[x,2])+' '+str(conn_matrix_char[x,3])+' '+str(conn_matrix_char[x,4])+' '+str(conn_matrix_char[x,5])+'\n')
     file.close()

  ###Now I have to give all the atoms the correct labeling based on the input script
     oneprint(myid,'Typing atoms...')
     #print(atypedict)
     ATOM_TYPE_REAL,LAMMPS_ATOMINDEXES,LAMMPS_CHARGES,conn_matrix_char_TYPED=types_charges(myid,natoms,atom_type,conn_matrix_nums,conn_matrix_char)
     #print(ATOM_TYPE_REAL)
     #print(ATOM_TYPE_REAL[0])
     #print(LAMMPS_ATOMINDEXES.dtype)
     #print(LAMMPS_CHARGES.dtype)
     #print(conn_matrix_char_TYPED.dtype)
     oneprint(myid,'Typed atoms!')

  ###Find all the number of organic linkers (center of masses,lists of atoms in each/nums&chars)
     oneprint(myid,'Linker connectivity...')
     NEGLECT_METALS_LIST=metals_indexes(myid,natoms,ATOM_TYPE_REAL,neglect_atype)
     print((NEGLECT_METALS_LIST))
     NCL,nclustersizes,ptr=percolate(myid,neglect_atype,ATOM_TYPE_REAL,conn_matrix_nums,NEGLECT_METALS_LIST)
     print('#LIGANDS::'+str(NCL).replace('[[','').replace(']]',''))
     #print(nclustersizes)
     #print(ptr)
     cluster_dict=linker_dict_creator(myid,natoms,nclustersizes,natomsinlinker,ptr)
     oneprint(myid,'Linker connectivity determined!')

     #Get the organic linker fragment 
     oneprint(myid,'Reading linker fragment.')
     file=open(input+'/'+partname,'r')
     atype_frag=[]
     xpos_frag=[]
     ypos_frag=[]
     zpos_frag=[]
     flag_frag=[]
     lines=file.readlines()
     for line in lines:
         #print(line)
         elements=line.split()
         atype_frag.append(elements[0])
         xpos_frag.append(elements[1])
         ypos_frag.append(elements[2])
         zpos_frag.append(elements[3])
         flag_frag.append(elements[4])
     file.close()
     xpos_frag=np.array(xpos_frag,dtype='f')
     ypos_frag=np.array(ypos_frag,dtype='f')
     zpos_frag=np.array(zpos_frag,dtype='f')
     PRIM_POS=flag_frag.index('PRIMARY')
     SECO_POS=flag_frag.index('SECONDARY')
     TERT_POS=flag_frag.index('TERTIARY')
     #print(PRIM_POS,SECO_POS,TERT_POS)

  ###Perform Reverse MC on a grid to decide which linkers to switch (base this on a short range order parameter) 
     oneprint(myid,'Deciding which linkers to switch...')

     #Set up the arrays for the reverse MC
     #Scan through the clusters to determine what is connected to what:

     #(1)do this by looking at the 2 Ns of a particlar cluster
     #(2)set the two Zn indexes it is connected to     
     METAL_DICT=defaultdict(list)
     for clust in range(NCL):
          ATOMLIST_TEST=cluster_dict[clust+1]
          Nflag=0
          for atom in range(natomsinlinker[0]):
              atomid=int(ATOMLIST_TEST[atom])
              if str(ATOM_TYPE_REAL[atomid])==str(['N']) and Nflag==0:#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                   CONNROW_numsA=conn_matrix_nums[atomid,:]-1
                   #print(CONNROW_numsA)
                   #Check to see if any of the metal indexes are in the CONNROW_numsA and store in a linker dictionary
                   metal_attached_to_linkerA=[i for i in CONNROW_numsA if i in NEGLECT_METALS_LIST]
                   METAL_DICT[clust+1].append(metal_attached_to_linkerA[0])
                   Nflag=1
              elif str(ATOM_TYPE_REAL[atomid])==str(['N']) and Nflag==1:#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                   CONNROW_numsB=conn_matrix_nums[atomid,:]-1
                   #print(CONNROW_numsB)
                   #Check to see if any of the metal indexes are in the CONNROW_numsB and store in a linker dictionary
                   metal_attached_to_linkerB=[i for i in CONNROW_numsB if i in NEGLECT_METALS_LIST]
                   METAL_DICT[clust+1].append(metal_attached_to_linkerB[0])   
     print(METAL_DICT)

     #(3)starting with one cluster, find the 3 other clusters connected to 1 Zn and same for other Zn (this makes the NN)
     linker_NNlist=np.ones((NCL,6),dtype='int')
     for clusta in range(NCL):
         MET_LISTA=METAL_DICT[clusta+1]
         counter = 0
         for metalindex in range(2):
             met_inA=[int(MET_LISTA[metalindex])]
             for clustb in range(NCL):
                 MET_LISTB=METAL_DICT[clustb+1]
                 #print(str([i for i in met_inA if i in MET_LISTB]))
                 if str([i for i in met_inA if i in MET_LISTB])!='[]' and clusta != clustb:
                     #print(clusta,counter,clustb+1)
                     linker_NNlist[clusta,counter]=(clustb+1)
                     counter+=1
     #print(linker_NNlist)

     PERFECT_NUMBER2SWITCH=int(round(frac2switch*float(NCL),0))
     oneprint(myid,'Pick which linkers to switch.')
     linkersswitching=[]
     for clust in range(NCL):
         if len(linkersswitching)==PERFECT_NUMBER2SWITCH:
           break
         if np.random.ranf() < frac2switch+0.2:
             linkersswitching.append(clust+1)
     print('#CLUSTERS CHOOSEN::'+str(len(linkersswitching)))
     print('FRACTION OF NEW LINKER::'+str(float(len(linkersswitching))/float(NCL)))
     print('CLUSTERS CHOOSEN::'+str(linkersswitching))
     linker_composition=(float(len(linkersswitching))/float(NCL)) #THIS HAS TO BE ZIF-90(PD) since I LOOK AWAY FROM ZIF-8(GOLD) LINKERS 

     file=open(output+'/COMPOSITION.data','w')
     file.write(str(float(len(linkersswitching))/float(NCL))+'\n')
     file.write(str(len(linkersswitching))+'\n')
     file.write(str(PERFECT_NUMBER2SWITCH)+'\n')
     file.write(str(NCL[0][0])+'\n')
     file.close()

     #(4)for a certain composition, choose which linkers get a 0(old linker) or a 1(new linker) and update the binary NN list
     linker_binaryNNlist=np.zeros((NCL,7),dtype='int')
     for clust in range(NCL):
         clustID=clust+1
         if clustID in linkersswitching:
             linker_binaryNNlist[clust,0]=1
         for clustnum in range(NCL):
             clustIDB=clustnum+1
             for neighbor in range(6):
                 if clustIDB==linker_NNlist[clust,neighbor] and (clustIDB in linkersswitching):
                     linker_binaryNNlist[clust,neighbor+1]=1
     print('FIRST LINKER BINARY LIST')
     print(linkersswitching)
     print(linker_binaryNNlist)
     #print(np.sum(linker_binaryNNlist))


     #(5)perform the MC scheme to get a certain short range order parameter [alpha=(1-Probability of 0 given 1)/(composition of 1)]:: no 0s next to 1s (alpha==1) and all 0s next to 1s (alpha==-1)
     #(6)switch the linkers that have a 1 attached to their index
     oneprint(myid,'Run MC.')
     SOP_old=measure_SOP(linker_composition,linker_binaryNNlist,n_linkerNNs)
     SOP_new=SOP_old
     SOP_accepted=SOP_old
     error_SOP_old=np.abs(SOP_old-set_SOP)
     error_SOP_saved=[]

     #for i in range(10000):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
     while (error_SOP_old > RMC_errorTOL): #probably set a certain number of interations
         print(SOP_accepted,error_SOP_old)
         #--1--#Pick two random linkers (1 is old and 2 is new)
         linkerflag=0
         while (linkerflag==0): #so we don't pick the same linker
             rand_linker_1=np.random.random_integers(1,NCL)
             if (rand_linker_1 in linkersswitching): #one of these has to be an 8 and the other a 90
                 linkerflag=0
             else:
                 linkerflag=1

         linkerflag=0
         while (linkerflag==0): #so we don't pick the same linker
             rand_linker_2=np.random.random_integers(1,NCL)
             if (rand_linker_2 in linkersswitching): #one of these has to be an 8 and the other a 90
                 linkerflag=1
             else:
                 linkerflag=0
         #print(rand_linker_1,rand_linker_2)

         #--2--#Switch the two linkers we picked
         linker_binaryNNlist_new,linkersswitching_new=switch_linkers(rand_linker_1,rand_linker_2,NCL,linker_NNlist,linkersswitching,n_linkerNNs,linker_binaryNNlist)
         #linker_binaryNNlist_new,linkersswitching_new=switch_linkers(rand_linker_1,rand_linker_2,NCL,linker_NNlist,linkersswitching,n_linkerNNs)
         #print(linker_binaryNNlist_new)
         #print(linkersswitching_new)
           
         #--3--#Calculate the new SOP
         SOP_new=measure_SOP(linker_composition,linker_binaryNNlist_new,n_linkerNNs)

         #--4--#Run the MC acceptance criteria
         error_SOP_new=np.abs(SOP_new-set_SOP)
         if error_SOP_new < error_SOP_old:
             #accept the move
             #print('A')
             linker_binaryNNlist=linker_binaryNNlist_new
             error_SOP_old=error_SOP_new
             linkersswitching=linkersswitching_new
             error_SOP_saved.append(error_SOP_old)
             SOP_accepted=SOP_new
         elif error_SOP_new >= error_SOP_old:
             randaccept=np.random.ranf()
             if randaccept < np.exp(-MCbeta*error_SOP_new):
                 #accept the move
                 #print('A')
                 linker_binaryNNlist=linker_binaryNNlist_new
                 error_SOP_old=error_SOP_new
                 linkersswitching=linkersswitching_new
                 error_SOP_saved.append(error_SOP_old)
                 SOP_accepted=SOP_new
             elif randaccept >= np.exp(-MCbeta*error_SOP_new):
                 #reject the move
                 #print('R')
                 linker_binaryNNlist=linker_binaryNNlist
                 error_SOP_old=error_SOP_old
                 linkersswitching=linkersswitching
                 error_SOP_saved.append(error_SOP_old)
                 SOP_accepted=SOP_accepted
         #if error_SOP_old <= 0.02:
             #break
     print(SOP_new)
     print(error_SOP_old)
     #print(linkersswitching_new)
     #print(linker_binaryNNlist_new)
     #print(np.sum(linker_binaryNNlist))
     file=open(output+'/reverseMC.dat','w')
     for i in range(len(error_SOP_saved)):
         file.write(str(i+1)+' '+str(error_SOP_saved[i])+'\n')
     file.close()

     file=open(output+'/FINAL_SRO.data','w')
     file.write(str(SOP_new)+'\n')
     file.write(str(set_SOP)+'\n')
     file.write(str(error_SOP_old)+'\n')
     file.write(str(RMC_errorTOL)+'\n')
     file.close()

  ###Switch out certain linkers in the appropriate ratio
     oneprint(myid,'Swapping out linkers...')
     atomid_stored=[]
     #These did not have the right type until now
     xpos=np.array(xpos,dtype='f')
     ypos=np.array(ypos,dtype='f')
     zpos=np.array(zpos,dtype='f')
     for clust in range(len(linkersswitching)):
          CLUSTID=linkersswitching[clust]
          ATOMLIST_TEST=cluster_dict[CLUSTID]
          atype_org=[]
          xpos_org=[]
          ypos_org=[]
          zpos_org=[]

          xpos_CHECKING=[]
          ypos_CHECKING=[]
          zpos_CHECKING=[]

          #NEED TO FIND THE LOCATION OF THE PRIMARY CARBON ON THE ORIGINAL RING
          for atom in range(natomsinlinker[0]):
               atomid=int(ATOMLIST_TEST[atom])
               atomid_stored.append(atomid)
               atype_org.append(ATOM_TYPE_REAL[atomid,0])
              
               REALpos_vec=np.array((xpos[atomid],ypos[atomid],zpos[atomid]),dtype='f')
               REALpos_vec=REALpos_vec.reshape(1,3)

               if str(ATOM_TYPE_REAL[atomid])==str([primary_carbonORG]):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                   C1x_FIXED=float(REALpos_vec[0,0])
                   C1y_FIXED=float(REALpos_vec[0,1])
                   C1z_FIXED=float(REALpos_vec[0,2])

          Nflag=0
          for atom in range(natomsinlinker[0]):
               atomid=int(ATOMLIST_TEST[atom])
               atomid_stored.append(atomid)
               atype_org.append(ATOM_TYPE_REAL[atomid,0])
              
               REALpos_vec=np.array((xpos[atomid],ypos[atomid],zpos[atomid]),dtype='f')
               REALpos_vec=REALpos_vec.reshape(1,3)

               xpos_CHECKING.append(float(REALpos_vec[0,0]))
               ypos_CHECKING.append(float(REALpos_vec[0,1]))
               zpos_CHECKING.append(float(REALpos_vec[0,2]))

               REALpos_vec=np.array(REALpos_vec.T)
               #print(REALpos_vec.shape)
               FRACpos_vec=np.dot(C_mat_inv,REALpos_vec)
               #print(FRACpos_vec)
               #print(REALpos_vec)

               if np.sqrt((C1y_FIXED-float(REALpos_vec[1,0]))**2) > 7:
                   if float(REALpos_vec[1,0]) < float(C1y_FIXED):
                       y_add=y_expand
                       FRACpos_vec[1,0]=FRACpos_vec[1,0]+y_add
                   else:
                       y_add=y_expand
                       FRACpos_vec[1,0]=FRACpos_vec[1,0]-y_add

               elif np.sqrt((C1z_FIXED-float(REALpos_vec[2,0]))**2) > 7:
                   if float(REALpos_vec[2,0]) < float(C1z_FIXED):
                       z_add=z_expand
                       FRACpos_vec[2,0]=FRACpos_vec[2,0]+z_add
                   else:
                       z_add=z_expand
                       FRACpos_vec[2,0]=FRACpos_vec[2,0]-z_add

               elif np.sqrt((C1x_FIXED-float(REALpos_vec[0,0]))**2) > 7:
                   if float(REALpos_vec[0,0]) < float(C1x_FIXED):
                       x_add=x_expand
                       FRACpos_vec[0,0]=FRACpos_vec[0,0]+x_add
                   else:
                       x_add=x_expand
                       FRACpos_vec[0,0]=FRACpos_vec[0,0]-x_add

               #print(FRACpos_vec)
               REALpos_vec=np.dot(C_mat,FRACpos_vec)
               #print(REALpos_vec)

               xpos_org.append(float(REALpos_vec[0,0]))
               ypos_org.append(float(REALpos_vec[1,0]))
               zpos_org.append(float(REALpos_vec[2,0]))
               if str(ATOM_TYPE_REAL[atomid])==str([primary_carbonORG]):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                   C1x_or=float(REALpos_vec[0,0])
                   C1y_or=float(REALpos_vec[1,0])
                   C1z_or=float(REALpos_vec[2,0])
               if str(ATOM_TYPE_REAL[atomid])==str([primary_nitrogenORG1]) and Nflag==0:#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                   N1x_or=float(REALpos_vec[0,0])
                   N1y_or=float(REALpos_vec[1,0])
                   N1z_or=float(REALpos_vec[2,0])
                   Nflag=1
               if str(ATOM_TYPE_REAL[atomid])==str([primary_nitrogenORG2]) and Nflag==1:#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                   N2x_or=float(REALpos_vec[0,0])
                   N2y_or=float(REALpos_vec[1,0])
                   N2z_or=float(REALpos_vec[2,0])
          
          file=open(output+'/linkercheck_org.xyz','w')
          for i in range(natomsinlinker[0]):
              file.write('O '+str(xpos_org[i])+' '+str(ypos_org[i])+' '+str(zpos_org[i])+'\n')
          file.close()

          file=open(output+'/linkercheck_new.xyz','w')
          for i in range(natomsinlinker[0]):
              file.write('Cl '+str(xpos_CHECKING[i])+' '+str(ypos_CHECKING[i])+' '+str(zpos_CHECKING[i])+'\n')
          file.close()

          C1or=(C1x_or,C1y_or,C1z_or)
          N1or=(N1x_or,N1y_or,N1z_or)
          N2or=(N2x_or,N2y_or,N2z_or)
          NNor=midpoint(N1or,N2or)

          #if CLUSTID==10:
              #print(C1or,N1or,N2or)

          #ROTATE
          #___________________________________________________
          #FIRST ROTATION TO MATCH UP THE PLANES
          #___________________________________________________
          C1=(xpos_frag[PRIM_POS],ypos_frag[PRIM_POS],zpos_frag[PRIM_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          N1=(xpos_frag[SECO_POS],ypos_frag[SECO_POS],zpos_frag[SECO_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          N2=(xpos_frag[TERT_POS],ypos_frag[TERT_POS],zpos_frag[TERT_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          NN=midpoint(N1,N2)
          replaceplane_vec=plane2normal(C1,N1,N2) 
          origplane_vec=plane2normal(C1or,N1or,N2or)
          angle1=angle_between(replaceplane_vec,origplane_vec)
          axis=axisfinder(replaceplane_vec,origplane_vec)
          #print(angle1)
          xstore=[]
          ystore=[]
          zstore=[]
          for atom in range(len(atype_frag)):          
              x1,y1,z1=YPR_transform(xpos_frag[atom],ypos_frag[atom],zpos_frag[atom],(10,10,10),axis,np.degrees(angle1))
              xstore.append(x1)
              ystore.append(y1)
              zstore.append(z1)
              #atom_type=np.append(atom_type,atype_frag[atom])
              #xpos=np.append(xpos,x1)
              #ypos=np.append(ypos,y1)
              #zpos=np.append(zpos,z1)

          #TRANSLATE AT THE END OF ROTATION
          DISTx=C1x_or-float(xstore[PRIM_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          DISTy=C1y_or-float(ystore[PRIM_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          DISTz=C1z_or-float(zstore[PRIM_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          xstore=np.array(xstore,dtype='f')
          ystore=np.array(ystore,dtype='f')
          zstore=np.array(zstore,dtype='f')
          xstore=xstore+DISTx
          ystore=ystore+DISTy
          zstore=zstore+DISTz

          #___________________________________________________
          #SECOND ROTATION TO MATCH UP THE DIRECTIONALITY OF THE RINGS
          #___________________________________________________
          C1=(xstore[PRIM_POS],ystore[PRIM_POS],zstore[PRIM_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          N1=(xstore[SECO_POS],ystore[SECO_POS],zstore[SECO_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          N2=(xstore[TERT_POS],ystore[TERT_POS],zstore[TERT_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          NN=midpoint(N1,N2)
          #print(C1,N1,N2,NN)
          #print(C1or,N1or,N2or,NNor)
          #print(C1[0],N1[0])
          xC1=copy.deepcopy(float(C1[0]))
          yC1=copy.deepcopy(float(C1[1]))
          zC1=copy.deepcopy(float(C1[2]))
          xN1=copy.deepcopy(float(N1[0]))
          yN1=copy.deepcopy(float(N1[1]))
          zN1=copy.deepcopy(float(N1[2]))
          xC1or=copy.deepcopy(float(C1or[0]))
          yC1or=copy.deepcopy(float(C1or[1]))
          zC1or=copy.deepcopy(float(C1or[2]))
          xN1or=copy.deepcopy(float(N1or[0]))
          yN1or=copy.deepcopy(float(N1or[1]))
          zN1or=copy.deepcopy(float(N1or[2]))
          #print(xC-xN)
          NEWvec=(xC1-xN1,yC1-yN1,zC1-zN1)
          OLDvec=(xC1or-xN1or,yC1or-yN1or,zC1or-zN1or)
          #print(OLDvec,NEWvec)
          angle2=angle_between(NEWvec,OLDvec)
          axis=axisfinder(NEWvec,OLDvec)
          #print(np.degrees(angle2))
          xstore2=[]
          ystore2=[]
          zstore2=[]
          for atom in range(len(atype_frag)):          
              x1,y1,z1=YPR_transform(xstore[atom],ystore[atom],zstore[atom],(10,10,10),axis,np.degrees(angle2))
              xstore2.append(x1)
              ystore2.append(y1)
              zstore2.append(z1)
              #atom_type=np.append(atom_type,atype_frag[atom])
              #xpos=np.append(xpos,x1)
              #ypos=np.append(ypos,y1)
              #zpos=np.append(zpos,z1)

          #TRANSLATE AT THE END OF ROTATION
          DISTx=C1x_or-float(xstore2[PRIM_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          DISTy=C1y_or-float(ystore2[PRIM_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          DISTz=C1z_or-float(zstore2[PRIM_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
          xstore2=np.array(xstore2,dtype='f')
          ystore2=np.array(ystore2,dtype='f')
          zstore2=np.array(zstore2,dtype='f')
          xstore2=xstore2+DISTx
          ystore2=ystore2+DISTy
          zstore2=zstore2+DISTz

          if 'NO' == 'YES':
              #___________________________________________________
              #THIRD ROTATION TO TILT THE LARGE RINGS SO NO OVERLAPPING ATOMS
              #___________________________________________________
              C1=(xstore2[PRIM_POS],ystore2[PRIM_POS],zstore2[PRIM_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
              N1=(xstore2[SECO_POS],ystore2[SECO_POS],zstore2[SECO_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
              N2=(xstore2[TERT_POS],ystore2[TERT_POS],zstore2[TERT_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
              #print(N1,N2)
              xN1=copy.deepcopy(float(N1[0]))
              yN1=copy.deepcopy(float(N1[1]))
              zN1=copy.deepcopy(float(N1[2]))
              xN2=copy.deepcopy(float(N2[0]))
              yN2=copy.deepcopy(float(N2[1]))
              zN2=copy.deepcopy(float(N2[2]))
              angle2=-20
              xstore3TEMPA=[]
              ystore3TEMPA=[]
              zstore3TEMPA=[]
              for atom in range(len(atype_frag)):          
                  x1,y1,z1=YPR_transform2(xstore2[atom],ystore2[atom],zstore2[atom],(xN1,yN1,zN1),(xN2,yN2,zN2),angle2)
                  xstore3TEMPA.append(x1)
                  ystore3TEMPA.append(y1)
                  zstore3TEMPA.append(z1)

              C1=(xstore2[PRIM_POS],ystore2[PRIM_POS],zstore2[PRIM_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
              N1=(xstore2[SECO_POS],ystore2[SECO_POS],zstore2[SECO_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
              N2=(xstore2[TERT_POS],ystore2[TERT_POS],zstore2[TERT_POS])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
              #print(N1,N2)
              xN1=copy.deepcopy(float(N1[0]))
              yN1=copy.deepcopy(float(N1[1]))
              zN1=copy.deepcopy(float(N1[2]))
              xN2=copy.deepcopy(float(N2[0]))
              yN2=copy.deepcopy(float(N2[1]))
              zN2=copy.deepcopy(float(N2[2]))
              angle2=20
              xstore3TEMPB=[]
              ystore3TEMPB=[]
              zstore3TEMPB=[]
              for atom in range(len(atype_frag)):          
                  x1,y1,z1=YPR_transform2(xstore2[atom],ystore2[atom],zstore2[atom],(xN1,yN1,zN1),(xN2,yN2,zN2),angle2)
                  xstore3TEMPB.append(x1)
                  ystore3TEMPB.append(y1)
                  zstore3TEMPB.append(z1)

              if (xstore3TEMPB[1]-8.5)**2+(ystore3TEMPB[1]-8.5)**2+(zstore3TEMPB[1]-8.5)**2 < (xstore3TEMPA[1]-8.5)**2+(ystore3TEMPA[1]-8.5)**2++(zstore3TEMPA[1]-8.5)**2:
                  xstore3=xstore3TEMPA
                  ystore3=ystore3TEMPA
                  zstore3=zstore3TEMPA
              else:
                  xstore3=xstore3TEMPB
                  ystore3=ystore3TEMPB
                  zstore3=zstore3TEMPB


          #for atom in range(len(atype_org)):
              #atom_type=np.append(atom_type,atype_org[atom])
              #xpos=np.append(xpos,xpos_org[atom])
              #ypos=np.append(ypos,ypos_org[atom])
              #zpos=np.append(zpos,zpos_org[atom])

          for atom in range(len(atype_frag)):
              atom_type=np.append(atom_type,atype_frag[atom])
              xpos=np.append(xpos,xstore2[atom])
              ypos=np.append(ypos,ystore2[atom])
              zpos=np.append(zpos,zstore2[atom])

          #for atom in range(len(atype_frag)):
              #atom_type=np.append(atom_type,atype_frag[atom])
              #xpos=np.append(xpos,xstore3[atom])
              #ypos=np.append(ypos,ystore3[atom])
              #zpos=np.append(zpos,zstore3[atom])

     #Remove that row from these arrays
     atom_type=np.delete(atom_type,atomid_stored)
     xpos=np.delete(xpos,atomid_stored)
     ypos=np.delete(ypos,atomid_stored)
     zpos=np.delete(zpos,atomid_stored)
#------------------------------------------------------------
     file=open(output+'/'+str(name)+'TESTREMOVE.xyz','w')
     file.write(str(len(atom_type))+'\n')
     file.write(str(name)+' '+str(aLEN)+' '+str(bLEN)+' '+str(cLEN)+' '+str(alphaD)+' '+str(betaD)+' '+str(gammaD)+'\n')
     natomsNOW=len(atom_type)
     storeCOORDSFORFUN=np.zeros([natomsNOW,3],dtype='f')
     for i in range(len(atom_type)):
         REALpos_vec=np.array((xpos[i],ypos[i],zpos[i]),dtype='f')
         REALpos_vec=REALpos_vec.reshape(1,3)
         REALpos_vec=np.array(REALpos_vec.T)
         #print(REALpos_vec.shape)
         FRACpos_vec=np.dot(C_mat_inv,REALpos_vec)
         #print(FRACpos_vec)

         #THIS MAY NEED TO BE FIXED
         if float(FRACpos_vec[0,0]) > float(x_expand):
             x_add=-(x_expand)
             FRACpos_vec[0,0]=FRACpos_vec[0,0]+x_add
         elif float(FRACpos_vec[0,0]) < float(0):
             x_add=x_expand
             FRACpos_vec[0,0]=FRACpos_vec[0,0]+x_add
         else:
             x_add=0
         if float(FRACpos_vec[1,0]) > float(y_expand):
             y_add=-(y_expand)
             FRACpos_vec[1,0]=FRACpos_vec[1,0]+y_add
         elif float(FRACpos_vec[1,0]) < float(0):
             y_add=y_expand
             FRACpos_vec[1,0]=FRACpos_vec[1,0]+y_add
         else:
             y_add=0
         if float(FRACpos_vec[2,0]) > float(z_expand):
             z_add=-(z_expand)
             FRACpos_vec[2,0]=FRACpos_vec[2,0]+z_add
         elif float(FRACpos_vec[2,0]) < float(0):
             z_add=z_expand
             FRACpos_vec[2,0]=FRACpos_vec[2,0]+z_add
         else:
             z_add=0
         REALpos_vec=np.dot(C_mat,FRACpos_vec)
         storeCOORDSFORFUN[i,0]=FRACpos_vec[0,0]
         storeCOORDSFORFUN[i,1]=FRACpos_vec[1,0]
         storeCOORDSFORFUN[i,2]=FRACpos_vec[2,0]

         file.write(str(atom_type[i]).replace("['",'').replace("']",'')+' '+str(float(REALpos_vec[0,0]))+' '+str(float(REALpos_vec[1,0]))+' '+str(float(REALpos_vec[2,0]))+'\n')
     file.write('\r')
     file.close()

     returnindex8888=cif_write(myid,'testremove',natomsNOW,atom_type,storeCOORDSFORFUN,aLEN,bLEN,cLEN,alphaD,betaD,gammaD,x_expand,y_expand,z_expand)

     coord_reg_HYBRID=np.vstack([xpos,ypos,zpos]) 
     coord_reg_HYBRID=coord_reg_HYBRID.T
     frac_HYBRIDLARGE=np.dot(C_mat_inv_LARGE,coord_reg_HYBRID.T)
     cif_writeALTERNATE(myid,'correctexpandedhybrid',natomsNOW,atom_type,frac_HYBRIDLARGE,aLEN,bLEN,cLEN,alphaD,betaD,gammaD,x_expand,y_expand,z_expand)

     #time.sleep(100)
#------------------------------------------------------------
  ###REGENERATE THE CONNECTIVITY
     if 'YES'=='YES':
         oneprint(myid,'Regenerating the conectivity.')
         #(7)regenerate the connectivity, typing of all the atoms, and get the new linker dictionary
         #Get fractional coordinates for the hybrid ZIF
         coord_reg_HYBRID=np.vstack([xpos,ypos,zpos]) 
         coord_reg_HYBRID=coord_reg_HYBRID.T
         coord_frac_HYBRID=np.dot(C_mat_inv,coord_reg_HYBRID.T)
         coord_frac_HYBRID=np.array(coord_frac_HYBRID.T)
         natoms=int(len(atom_type))
         print('NATOMS:::::::',natoms)
         for pos in range(len(atom_type)):#---------------------------------------------------------------------------------------------------------------------------> MAKE A FUNCTION FOR READABILITY
             if coord_frac_HYBRID[pos,0] > x_expand:
                 coord_frac_HYBRID[pos,0]=coord_frac_HYBRID[pos,0]-x_expand
             elif coord_frac_HYBRID[pos,0] < 0:
                 coord_frac_HYBRID[pos,0]=coord_frac_HYBRID[pos,0]+x_expand

             if coord_frac_HYBRID[pos,1] > y_expand:
                 coord_frac_HYBRID[pos,1]=coord_frac_HYBRID[pos,1]-y_expand
             elif coord_frac_HYBRID[pos,1] < 0:
                 coord_frac_HYBRID[pos,1]=coord_frac_HYBRID[pos,1]+y_expand

             if coord_frac_HYBRID[pos,2] > z_expand:
                 coord_frac_HYBRID[pos,2]=coord_frac_HYBRID[pos,2]-z_expand
             elif coord_frac_HYBRID[pos,2] < 0:
                 coord_frac_HYBRID[pos,2]=coord_frac_HYBRID[pos,2]+z_expand

         #print(coord_frac_HYBRID.shape)
         coord_reg_HYBRID=np.dot(C_mat,coord_frac_HYBRID.T)
         coord_reg_HYBRID=np.array(coord_reg_HYBRID.T)
         #print(coord_reg_HYBRID.shape)

         atom_type=np.array(atom_type,dtype='S7')
         natoms=len(atom_type)
         atom_typeINT=[]
         for atom in range(natoms):
             if atom_type[atom]=='H':
                 atom_typeINT.append(1)
             if atom_type[atom]=='C':
                 atom_typeINT.append(2)
             if atom_type[atom]=='N':
                 atom_typeINT.append(3)
             if atom_type[atom]=='O':
                 atom_typeINT.append(4)
             if atom_type[atom]=='Zn':
                 atom_typeINT.append(5)
         atom_typeINT=np.array(atom_typeINT,dtype='int')

         conn_matrix1=np.zeros([natoms,13])
         conn_matrix_char=np.zeros([natoms,6],dtype='a3')
         conn_matrix_nums1=np.zeros([natoms,6],dtype='int')
         bond_fix = 0.4
         radii_array=[]
         for atom in range(natoms):
             atom1=atom_type[atom]
             rad1=cov_radius_assignment(atype,radii,atom1)
             radii_array.append(rad1)
         radii_array1=np.array(radii_array,dtype='f')
         #print (CONNECT.connectivity_mapping.__doc__)
         #print('sleeping')
         #time.sleep(30)
         
         conn_matrix1,conn_matrix_nums1=COMBINEDv2.connectivity_mapping(atom_typeINT,radii_array1,coord_frac_HYBRID,C_mat,x_expand,y_expand,z_expand,conn_matrix1,conn_matrix_nums1,bond_corr)
         #conn_matrix,conn_matrix_nums,conn_matrix_char=connectivity_mapping(natoms,radii,atom_type,atype,coord_frac,C_mat,x_expand,y_expand,z_expand) 
         #print(conn_matrix_char)
         #print(conn_matrix_nums1)
         #print(conn_matrix1)
         conn_matrix=copy.deepcopy(conn_matrix1)
         conn_matrix_nums=copy.deepcopy(conn_matrix_nums1)
         for atom in range(natoms):
             for i in range(6):
                 if conn_matrix_nums[atom,i] != 0:
                     conn_matrix_char[atom,i]=atom_type[int(conn_matrix_nums[atom,i])-1]
         #print(conn_matrix_char)

         file=open(output+'/connect_check3_char.xyz','w')
         for x in range(natoms):
             file.write(str(conn_matrix_char[x,0])+' '+str(conn_matrix_char[x,1])+' '+str(conn_matrix_char[x,2])+' '+str(conn_matrix_char[x,3])+' '+str(conn_matrix_char[x,4])+' '+str(conn_matrix_char[x,5])+'\n')
         file.close()

         file=open(output+'/connect_check4_nums.xyz','w')
         for x in range(natoms):
             file.write(str(conn_matrix_nums[x,0])+' '+str(conn_matrix_nums[x,1])+' '+str(conn_matrix_nums[x,2])+' '+str(conn_matrix_nums[x,3])+' '+str(conn_matrix_nums[x,4])+' '+str(conn_matrix_nums[x,5])+'\n')
         file.close()

  ###RETYPE ALL THE ATOMS 
         #print(atypedict)
         ATOM_TYPE_REAL_HYBRID,LAMMPS_ATOMINDEXES_HYBRID,LAMMPS_CHARGES_HYBRID,conn_matrix_char_TYPED_HYBRID=types_charges(myid,natoms,atom_type,conn_matrix_nums,conn_matrix_char)
         print(conn_matrix_char_TYPED_HYBRID)

         #coord_frac_HYBRIDLARGE=(np.dot(C_mat_inv_LARGE,coord_reg_HYBRID.T))
         #print(coord_frac_HYBRIDLARGE.shape)
         #coord_frac_HYBRIDLARGE=np.array(coord_frac_HYBRIDLARGE.T)
         #returnindex8888=cif_write(myid,'correctexpandedhybrid_TYPED',natoms,ATOM_TYPE_REAL_HYBRID,coord_frac_HYBRIDLARGE,aLEN,bLEN,cLEN,alphaD,betaD,gammaD,x_expand,y_expand,z_expand)
         returnindex8888=cif_writeALTERNATE(myid,'correctexpandedhybrid_TYPED',natoms,ATOM_TYPE_REAL_HYBRID,frac_HYBRIDLARGE,aLEN,bLEN,cLEN,alphaD,betaD,gammaD,x_expand,y_expand,z_expand)


  ###RE-PERCOLATE
         NEGLECT_METALS_LIST_HYBRID=metals_indexes(myid,natoms,ATOM_TYPE_REAL_HYBRID,neglect_atype)
         print('ZN to neglect:',len(NEGLECT_METALS_LIST_HYBRID))
         NCL_hyb,nclustersizes_hyb,ptr_hyb=percolate(myid,neglect_atype,ATOM_TYPE_REAL_HYBRID,conn_matrix_nums,NEGLECT_METALS_LIST_HYBRID)
         print('#LIGANDS::'+str(NCL_hyb).replace('[[','').replace(']]',''))
         print('CLUSTER SIZES',len(nclustersizes_hyb))
         #print((nclustersizes_hyb))
         #print(ptr_hyb)
         #time.sleep(100000)
         cluster_dict_hyb=linker_dict_creator(myid,natoms,nclustersizes_hyb,natomsinlinker,ptr_hyb)
         #print(cluster_dict_hyb)

         print('DONE PERCOLATING.')
         #print(cluster_dict_hyb)


  ### CLASSIFY THE WINDOWS IN THE STRUCTURE
     if 'YES'==intensity_flag:
         print('WINDOW HISTOGRAMING ALGORITHM STARTING')
         #print(cluster_dict_hyb) 

         METAL_DICT_HYBRID=defaultdict(list)
         for clust in range(NCL):
              ATOMLIST_TEST=cluster_dict_hyb[clust+1]
              Nflag=0
              for atom in range(len(ATOMLIST_TEST)):
                  atomid=int(ATOMLIST_TEST[atom])
                  if str(ATOM_TYPE_REAL_HYBRID[atomid])==str(['N']) and Nflag==0:#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                       CONNROW_numsA=conn_matrix_nums[atomid,:]-1
                       #print(CONNROW_numsA)
                       #Check to see if any of the metal indexes are in the CONNROW_numsA and store in a linker dictionary
                       metal_attached_to_linkerA=[i for i in CONNROW_numsA if i in NEGLECT_METALS_LIST_HYBRID]
                       METAL_DICT_HYBRID[clust+1].append(metal_attached_to_linkerA[0])
                       Nflag=1
                  elif str(ATOM_TYPE_REAL_HYBRID[atomid])==str(['N']) and Nflag==1:#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                       CONNROW_numsB=conn_matrix_nums[atomid,:]-1
                       #print(CONNROW_numsB)
                       #Check to see if any of the metal indexes are in the CONNROW_numsB and store in a linker dictionary
                       metal_attached_to_linkerB=[i for i in CONNROW_numsB if i in NEGLECT_METALS_LIST_HYBRID]
                       METAL_DICT_HYBRID[clust+1].append(metal_attached_to_linkerB[0])   
         #print(METAL_DICT_HYBRID) #this shows all the ZN atoms that are attached by linkers

         file=open(output+'/ZNs.xyz','w')
         for Zn in NEGLECT_METALS_LIST_HYBRID:
           file.write('Zn'+str(Zn)+' '+str(coord_reg_HYBRID[Zn,0])+' '+str(coord_reg_HYBRID[Zn,1])+' '+str(coord_reg_HYBRID[Zn,2])+'\n')
         file.close()

         linker_NNlist=np.ones((NCL,6),dtype='int')
         for clusta in range(NCL):
             MET_LISTA=METAL_DICT_HYBRID[clusta+1]
             counter = 0
             for metalindex in range(2):
                 met_inA=[int(MET_LISTA[metalindex])]
                 for clustb in range(NCL):
                     MET_LISTB=METAL_DICT_HYBRID[clustb+1]
                         #print(str([i for i in met_inA if i in MET_LISTB]))
                     if str([i for i in met_inA if i in MET_LISTB])!='[]' and clusta != clustb:
                         #print(clusta,counter,clustb+1)
                         linker_NNlist[clusta,counter]=(clustb+1)
                         counter+=1
         #print(linker_NNlist)

         COM_cluster_TYPES=defaultdict(list)
         #CARTESIAN COORDINATES
         for clust in range(NCL_hyb):
             clust_atom_list=cluster_dict_hyb[clust+1]
             linkerflag=0
             for atom in range(len(clust_atom_list)):
                 atomid=int(clust_atom_list[atom])
                 if str(ATOM_TYPE_REAL_HYBRID[atomid])==str(['H3']):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                     linkerflag=0
                     COM_cluster_TYPES[clust+1].append(linkerflag)
                     break
                 if str(ATOM_TYPE_REAL_HYBRID[atomid])==str(['H4']):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                     linkerflag=1
                     COM_cluster_TYPES[clust+1].append(linkerflag)
                     break
         #print(COM_cluster_TYPES)
           

         molecule,ZnNNdict=linkerconnection_format(METAL_DICT_HYBRID)
         print('CONNECTIVITY: ',molecule)
         #ring=getring(molecule, 40)
         rings=getrings(molecule)
         print('NUMBER OF RINGS: ',len(rings))
         #time.sleep(100)
         #print(rings)
         RINGTYPES=postprocessrings(rings,METAL_DICT_HYBRID,ZnNNdict,COM_cluster_TYPES)
         HISTOGRAMBIN,HISTOGRAMCOUNTS=histo(rings,RINGTYPES)
         file=open(output+'/HISTOGRAM_LINKER_COUNTS.data','w')
         for h in range(len(HISTOGRAMCOUNTS)):
           file.write(str(HISTOGRAMCOUNTS[h])+'\n')
         file.close()
         file=open(output+'/HISTOGRAM_LINKERS.data','w')
         file.write('nsamples: '+str(len(HISTOGRAMCOUNTS))+'\n\n')
         for h in range(len(HISTOGRAMBIN)):
           file.write(str(h)+' '+str(HISTOGRAMBIN[h])+'\n')
         file.close()
         print('Histogram file written.')

         #time.sleep(100000)

  ###INTENSITY CURVES
         print('STARTING INTENSITY CURVE PROCEDURE.')
         #(8)calculate the COM of all the new linkers
         COM_cluster_dictreg=defaultdict(list)
         COM_cluster_dictfrac=defaultdict(list)
         #CARTESIAN COORDINATES
         clust_counter=0
         for clust in range(NCL_hyb):
             clust_atom_list=cluster_dict_hyb[clust+1]
             linkerflag=0
             for atom in range(len(clust_atom_list)):
                 atomid=int(clust_atom_list[atom])
                 if str(ATOM_TYPE_REAL_HYBRID[atomid])==str([MAJORITY_HYDROGEN]):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                     #RELOOP THROUGH ALL THE ATOMS OF THAT CLUSTER 
                     XPOS_HYDROGENS=[]
                     YPOS_HYDROGENS=[]
                     ZPOS_HYDROGENS=[]
                     for atom888 in range(len(clust_atom_list)):
                       atomid=int(clust_atom_list[atom888])

                       if str(ATOM_TYPE_REAL_HYBRID[atomid])==str([MAJORITY_HYDROGEN]):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                         XPOS_HYDROGENS.append(float(coord_reg_HYBRID[atomid,0]))
                         YPOS_HYDROGENS.append(float(coord_reg_HYBRID[atomid,1]))
                         ZPOS_HYDROGENS.append(float(coord_reg_HYBRID[atomid,2]))
                         COM_cluster_dictreg[clust_counter+1].append(float(coord_reg_HYBRID[atomid,0]))
                         COM_cluster_dictreg[clust_counter+1].append(float(coord_reg_HYBRID[atomid,1]))
                         COM_cluster_dictreg[clust_counter+1].append(float(coord_reg_HYBRID[atomid,2]))
                         COM_cluster_dictreg[clust_counter+1].append(0)
                         clust_counter+=1
                     break

             for atom in range(len(clust_atom_list)):
                 atomid=int(clust_atom_list[atom])
                 if str(ATOM_TYPE_REAL_HYBRID[atomid])==str([MINORITY_HYDROGEN]):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                     #RELOOP THROUGH ALL THE ATOMS OF THAT CLUSTER 
                     XPOS_HYDROGENS=[]
                     YPOS_HYDROGENS=[]
                     ZPOS_HYDROGENS=[]
                     for atom888 in range(len(clust_atom_list)):
                       atomid=int(clust_atom_list[atom888])

                       if str(ATOM_TYPE_REAL_HYBRID[atomid])==str([MINORITY_HYDROGEN]):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                         XPOS_HYDROGENS.append(float(coord_reg_HYBRID[atomid,0]))
                         YPOS_HYDROGENS.append(float(coord_reg_HYBRID[atomid,1]))
                         ZPOS_HYDROGENS.append(float(coord_reg_HYBRID[atomid,2]))
                         COM_cluster_dictreg[clust_counter+1].append(float(coord_reg_HYBRID[atomid,0]))
                         COM_cluster_dictreg[clust_counter+1].append(float(coord_reg_HYBRID[atomid,1]))
                         COM_cluster_dictreg[clust_counter+1].append(float(coord_reg_HYBRID[atomid,2]))
                         COM_cluster_dictreg[clust_counter+1].append(1)
                         clust_counter+=1
                     break

         #print(COM_cluster_dictreg)

         file=open(output+'/LINKER_POSITIONS_Hs.xyz','w')
         file.write(str(len(COM_cluster_dictreg))+' '+'\n\n')
         for i in range(len(COM_cluster_dictreg)):
           dicreg=COM_cluster_dictreg[i+1]
           if dicreg[3] == 0:
             file.write('H3'+str(i+1)+' '+str(dicreg[0])+' '+str(dicreg[1])+' '+str(dicreg[2])+'\n')
           if dicreg[3] == 1:
             file.write('O4'+str(i+1)+' '+str(dicreg[0])+' '+str(dicreg[1])+' '+str(dicreg[2])+'\n')
         file.close()

         #FRACTIONAL COORDINATES
         clust_counter=0
         for clust in range(NCL_hyb):
             clust_atom_list=cluster_dict_hyb[clust+1]
             linkerflag=0
             for atom in range(len(clust_atom_list)):
                 atomid=int(clust_atom_list[atom])
                 if str(ATOM_TYPE_REAL_HYBRID[atomid])==str([MAJORITY_HYDROGEN]):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                     #RELOOP THROUGH ALL THE ATOMS OF THAT CLUSTER 
                     XPOS_HYDROGENS=[]
                     YPOS_HYDROGENS=[]
                     ZPOS_HYDROGENS=[]
                     for atom888 in range(len(clust_atom_list)):
                       atomid=int(clust_atom_list[atom888])

                       if str(ATOM_TYPE_REAL_HYBRID[atomid])==str([MAJORITY_HYDROGEN]):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                         XPOS_HYDROGENS.append(float(coord_frac_HYBRID[atomid,0]))
                         YPOS_HYDROGENS.append(float(coord_frac_HYBRID[atomid,1]))
                         ZPOS_HYDROGENS.append(float(coord_frac_HYBRID[atomid,2]))
                         COM_cluster_dictfrac[clust_counter+1].append(float(coord_frac_HYBRID[atomid,0]))
                         COM_cluster_dictfrac[clust_counter+1].append(float(coord_frac_HYBRID[atomid,1]))
                         COM_cluster_dictfrac[clust_counter+1].append(float(coord_frac_HYBRID[atomid,2]))
                         COM_cluster_dictfrac[clust_counter+1].append(0)
                         clust_counter+=1
                     break

             for atom in range(len(clust_atom_list)):
                 atomid=int(clust_atom_list[atom])
                 if str(ATOM_TYPE_REAL_HYBRID[atomid])==str([MINORITY_HYDROGEN]):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                     #RELOOP THROUGH ALL THE ATOMS OF THAT CLUSTER 
                     XPOS_HYDROGENS=[]
                     YPOS_HYDROGENS=[]
                     ZPOS_HYDROGENS=[]
                     for atom888 in range(len(clust_atom_list)):
                       atomid=int(clust_atom_list[atom888])

                       if str(ATOM_TYPE_REAL_HYBRID[atomid])==str([MINORITY_HYDROGEN]):#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
                         XPOS_HYDROGENS.append(float(coord_frac_HYBRID[atomid,0]))
                         YPOS_HYDROGENS.append(float(coord_frac_HYBRID[atomid,1]))
                         ZPOS_HYDROGENS.append(float(coord_frac_HYBRID[atomid,2]))
                         COM_cluster_dictfrac[clust_counter+1].append(float(coord_frac_HYBRID[atomid,0]))
                         COM_cluster_dictfrac[clust_counter+1].append(float(coord_frac_HYBRID[atomid,1]))
                         COM_cluster_dictfrac[clust_counter+1].append(float(coord_frac_HYBRID[atomid,2]))
                         COM_cluster_dictfrac[clust_counter+1].append(1)
                         clust_counter+=1
                     break


     #(9)build the intensity (#1s over #0s) vs distance curves 
         FRACCOORD_FOR_COM_EXPANSION_small=np.zeros((len(COM_cluster_dictfrac),4),dtype='float')
         for clust in range(len(COM_cluster_dictfrac)):
             coordlist_fracs=COM_cluster_dictfrac[clust+1]
             #print(coordlist_fracs)
             FRACCOORD_FOR_COM_EXPANSION_small[clust,0]=coordlist_fracs[0]
             FRACCOORD_FOR_COM_EXPANSION_small[clust,1]=coordlist_fracs[1]
             FRACCOORD_FOR_COM_EXPANSION_small[clust,2]=coordlist_fracs[2]
             FRACCOORD_FOR_COM_EXPANSION_small[clust,3]=coordlist_fracs[3]
         #print((FRACCOORD_FOR_COM_EXPANSION_small))

         FRACCOORD_large,LINKERDESCRIP_large=expand_LINKERPOS(len(COM_cluster_dictfrac),FRACCOORD_FOR_COM_EXPANSION_small,1,1,1)#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
         #print(len(FRACCOORD_large))
         #print(len(LINKERDESCRIP_large))
         FRACCOORD_large=np.array(FRACCOORD_large.T)
         print(FRACCOORD_large.shape)
         REGCOORD_large=np.dot(C_mat,FRACCOORD_large.T)
         REGCOORD_large=REGCOORD_large.T
         print(REGCOORD_large.shape)

         file=open(output+'/'+str(name)+'TEST_EXPAND_forintensity.xyz','w')
         file.write(str(len(REGCOORD_large))+'\n')
         file.write(str(name)+' '+str(aLEN)+' '+str(bLEN)+' '+str(cLEN)+' '+str(alphaD)+' '+str(betaD)+' '+str(gammaD)+'\n')
         for i in range(len(REGCOORD_large)):
             if int(LINKERDESCRIP_large[i]) == 0:
                 file.write(str('C').replace("['",'').replace("']",'')+' '+str(float(REGCOORD_large[i,0]))+' '+str(float(REGCOORD_large[i,1]))+' '+str(float(REGCOORD_large[i,2]))+'\n')
             if int(LINKERDESCRIP_large[i]) == 1:
                 file.write(str('O').replace("['",'').replace("']",'')+' '+str(float(REGCOORD_large[i,0]))+' '+str(float(REGCOORD_large[i,1]))+' '+str(float(REGCOORD_large[i,2]))+'\n')
         file.write('\r')
         file.close()

         #FIND ALL THE 0 LINKERS
         allZEROclusters=[]
         for clust in range(len(COM_cluster_dictfrac)):
             coordlist_fracs=COM_cluster_dictfrac[clust+1]
             if int(coordlist_fracs[3])==0:
                 allZEROclusters.append(clust+1)

         print(allZEROclusters)
         #time.sleep(1000) 

#CHANGE THIS SECTION OF CODE##########################
     if 'YES'=='NO':
           #distances_to_look_from=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,30]#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
           #distances_to_look_from=np.arange(start=1,stop=34,step=0.2)#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
           #allZEROclusters=[1,2,3]
           reshapefactor=len(distances_to_look_from)
           intense_storedA=np.zeros((reshapefactor,len(allZEROclusters)),dtype='float')
           intense_storedB=np.zeros((reshapefactor,len(allZEROclusters)),dtype='float')

           for LINKER in range(len(allZEROclusters)):
               #STANDING ON ZIF-8 LINKERS
               stand_on_linker=int(allZEROclusters[LINKER])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
               STAND_COORD=COM_cluster_dictreg[stand_on_linker]
               print(STAND_COORD)
               xSTAND_COORD=float(STAND_COORD[0])
               ySTAND_COORD=float(STAND_COORD[1])
               zSTAND_COORD=float(STAND_COORD[2])

               for D in range(len(distances_to_look_from)):
                   DSQcheck=(float(distances_to_look_from[D]))**2
                   countLA=0
                   countLB=0
                   for clust in range(len(LINKERDESCRIP_large)):
                       rx=xSTAND_COORD-REGCOORD_large[clust,0]
                       rx=rx-(aLEN*x_expand)*round(rx/float(aLEN*x_expand))

                       ry=ySTAND_COORD-REGCOORD_large[clust,1]
                       ry=ry-(bLEN*y_expand)*round(ry/float(bLEN*y_expand))

                       rz=zSTAND_COORD-REGCOORD_large[clust,2]
                       rz=rz-(cLEN*z_expand)*round(rz/float(cLEN*z_expand))

                       Dsq=float(rx)**2+float(ry)**2+float(rz)**2

                       if Dsq <= DSQcheck:
                           if int(LINKERDESCRIP_large[clust])==1:
                               countLA+=1
                           elif int(LINKERDESCRIP_large[clust])==0:
                               countLB+=1
                   intense_storedA[D,LINKER]=countLA
                   intense_storedB[D,LINKER]=countLB

           #print(intense_storedA)
           LAAVG=intense_storedA.mean(axis=1)
           LBAVG=intense_storedB.mean(axis=1)
           intensity=LAAVG/(LAAVG+3*LBAVG)
           #print(intensity.shape)
           #print(LAAVG.shape)
           #print(LBAVG.shape)
           file=open(output+'/INTENSITY1.data','w')
           for index in range(len(LAAVG)):
               file.write(str(LAAVG[index])+' '+str(LBAVG[index])+' '+str(distances_to_look_from[index])+' '+str(intensity[index])+'\n')
           file.close()
#####################################
     if 'YES'==intensity_flag:
         #DISTANCES_stored=np.zeros((len(allZEROclusters),len(COM_cluster_dictfrac)-1),dtype='float') #ROWS ZIF8 #COLUMNS ALL OTHER LINKERS (DO NOT DOUBLE COUNT)
         #LINKERTYPES_stored=np.zeros((len(allZEROclusters),len(COM_cluster_dictfrac)-1),dtype='int') #ROWS ZIF8 #COLUMNS ALL OTHER LINKERS (DO NOT DOUBLE COUNT)

         distances_to_look_from=copy.deepcopy(distances_to_look_fromD)
         DISTANCES_stored=np.zeros((len(allZEROclusters),len(COM_cluster_dictfrac)),dtype='float') #ROWS ZIF8 #COLUMNS ALL OTHER LINKERS (DO NOT DOUBLE COUNT)
         LINKERTYPES_stored=np.zeros((len(allZEROclusters),len(COM_cluster_dictfrac)),dtype='int') #ROWS ZIF8 #COLUMNS ALL OTHER LINKERS (DO NOT DOUBLE COUNT)

         for LINKER in range(len(allZEROclusters)):
             #STANDING ON ZIF-8 LINKERS
             stand_on_linker=int(allZEROclusters[LINKER])#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
             STAND_COORD=COM_cluster_dictreg[stand_on_linker]
             #print(STAND_COORD)
             xSTAND_COORD=float(STAND_COORD[0])
             ySTAND_COORD=float(STAND_COORD[1])
             zSTAND_COORD=float(STAND_COORD[2])
             fakecounter=0
             for clust in range(len(LINKERDESCRIP_large)):
                     rx=xSTAND_COORD-REGCOORD_large[clust,0]
                     rx=rx-(aLEN*x_expand)*round(rx/float(aLEN*x_expand))

                     ry=ySTAND_COORD-REGCOORD_large[clust,1]
                     ry=ry-(bLEN*y_expand)*round(ry/float(bLEN*y_expand))

                     rz=zSTAND_COORD-REGCOORD_large[clust,2]
                     rz=rz-(cLEN*z_expand)*round(rz/float(cLEN*z_expand))

                     Dsq=float(rx)**2+float(ry)**2+float(rz)**2
                     #print(int(clust+1),int(stand_on_linker))
                     #if int(clust+1) != int(stand_on_linker):
                     if 'YES' == 'YES':
                       DISTANCES_stored[LINKER,fakecounter]=np.sqrt(Dsq)
                       if int(LINKERDESCRIP_large[clust])==1:
                         LINKERTYPES_stored[LINKER,fakecounter]=1
                       fakecounter+=1

         #print(DISTANCES_stored)
         #print(LINKERTYPES_stored)

         #distances_to_look_from=np.arange(start=1,stop=100,step=1)#---------------------------------------------------------------------------------------------------------------------------> DO IN INPUT SCRIPT
         #INTENSITYSTORE=np.zeros((len(distances_to_look_from),len(allZEROclusters),len(COM_cluster_dictfrac)-1),dtype='float')
         #PARTIALINTENSITYSTORE=np.zeros((len(distances_to_look_from),len(allZEROclusters),len(COM_cluster_dictfrac)-1),dtype='float')

         INTENSITYSTORE=np.zeros((len(distances_to_look_from),len(allZEROclusters),len(COM_cluster_dictfrac)),dtype='float')
         PARTIALINTENSITYSTORE=np.zeros((len(distances_to_look_from),len(allZEROclusters),len(COM_cluster_dictfrac)),dtype='float')

         for d in range(len(distances_to_look_from)):
           print(d)
           dist=distances_to_look_from[d]
           DSQcheck=dist**2
           for L in range(len(allZEROclusters)):
              for LINKER in range(len(LINKERDESCRIP_large)-1):
                 Dsq=DISTANCES_stored[L,LINKER]**2
                 #MAYBE I DON"T COUNT LINKERS NOT WITHIN THE SEARCH RADIUS (DOES NOT MAKE MUCH SENSE THOUGH...)
                 if Dsq <= DSQcheck: 
                   I=math.erfc(DISTANCES_stored[L,LINKER]/(dist*np.sqrt(math.pi)))
                   #I=1
                 else:
                   I=0
                 #I=math.erfc(DISTANCES_stored[L,LINKER]/(dist*math.pi))

                 if LINKERTYPES_stored[L,LINKER] == 1:
                   INTENSITYSTORE[d,L,LINKER]=I
                   PARTIALINTENSITYSTORE[d,L,LINKER]=I
                 elif LINKERTYPES_stored[L,LINKER] == 0:
                   INTENSITYSTORE[d,L,LINKER]=I
                   PARTIALINTENSITYSTORE[d,L,LINKER]=0

         #print(PARTIALINTENSITYSTORE)
         #print(INTENSITYSTORE)

         print('writing file finally')
         ROSS_DATA=[]
         file=open(output+'/INTENSITY2.data','w')
         for d in range(len(distances_to_look_from)):
            print(d)
            dist=distances_to_look_from[d]
            PARTIAL_AVG=np.sum(PARTIALINTENSITYSTORE[d,:,:])/float(len(allZEROclusters))
            TOTAL_AVG=np.sum(INTENSITYSTORE[d,:,:])/float(len(allZEROclusters))
            if float(TOTAL_AVG) > 0:
              INTENSITY_AVG=float(PARTIAL_AVG)/float(TOTAL_AVG)
            else:
              INTENSITY_AVG=0
            file.write(str(dist)+' '+str(INTENSITY_AVG)+' '+str(PARTIAL_AVG)+' '+str(TOTAL_AVG)+'\n')
            ROSS_DATA.append(INTENSITY_AVG)
         file.close()

         if DISTANCE_SPACING_TYPE=='USER_DEFINED':
           DEVabs=np.abs((KRISHNA_DATA-ROSS_DATA)/KRISHNA_DATA) #PROBLEM WITH THIS IS THAT THE BIG DEVIATIONS WILL DOMINATE THE SMALL ONES
           print(DEVabs)
           sumDEVabs=np.sum(DEVabs)/len(DEVabs)
           print('MAPE:',sumDEVabs*100)
           file=open(output+'/FINAL_MAPE.data','w')
           file.write(str(sumDEVabs*100)+'\n')
           file.close()

           file=open(output+'/PARITY_intensity.data','w')
           for i in range(len(ROSS_DATA)):
             file.write(str(ROSS_DATA[i])+' '+str(KRISHNA_DATA[i])+'\n')
           file.close()

         #time.sleep(100)

  ###Find all BONDS AND ANGLES
     if 'YES'=='YES':
         oneprint(myid,'Found bonds and angles.')
         bond_count_HYBRID,angle_count_HYBRID,BOND_ARRAY_HYBRID,ANGLE_ARRAY_HYBRID=bonds_angles(myid,natoms,conn_matrix_char,conn_matrix_nums)
         oneprint(myid,'#BONDS::'+str(bond_count_HYBRID))
         oneprint(myid,'#ANGLES::'+str(angle_count_HYBRID))
         #print(ANGLE_ARRAY)
         #print(BOND_ARRAY.dtype)
         #print(ANGLE_ARRAY.dtype)

  ###Find all PROPER DIHEDRALS
         dihedral_count_HYBRID,DIHEDRAL_ARRAY_HYBRID=dihedrals(myid,natoms,conn_matrix_char,conn_matrix_nums,BOND_ARRAY_HYBRID)
         oneprint(myid,'#DIHEDRALS::'+str(dihedral_count_HYBRID))


  ###WRITE LAMMPS DATA FILE (regular coordinates)
     if 'YES'=='YES':
         #This function includes keeping only certain dihedrals
         print('WRITING LAMMPS DATA FILE')
         lammps_write(myid,name,natoms,ATOM_TYPE_REAL_HYBRID,LAMMPS_ATOMINDEXES_HYBRID,LAMMPS_CHARGES_HYBRID,coord_reg_HYBRID,BOND_ARRAY_HYBRID,ANGLE_ARRAY_HYBRID,DIHEDRAL_ARRAY_HYBRID,check_lmp_typing_flag,conn_matrix_char_TYPED_HYBRID,conn_matrix_nums)

  ###WRITE THE LJ PARAMETERS FILE
     if 'YES'=='YES':
         ljparameters(epsilon,sigma)

     print('YOU ARE DONE!!!')

  ###WRITE RASPA-2.0 CIF FILE (fractional coordinates)
     if 'YES'=='NO':
         returnindex1=cif_write(myid,name,natoms,ATOM_TYPE_REAL_HYBRID,coord_frac_HYBRID,aLEN,bLEN,cLEN,alphaD,betaD,gammaD,1,1,1)

         natoms_HYBRID=len(atom_type)
         name_alt=name+'_HYBRID'
         returnindex2=cif_write(myid,name_alt,natoms_HYBRID,atom_type,coord_frac_HYBRID,aLEN,bLEN,cLEN,alphaD,betaD,gammaD,1,1,1)

         HYBRIDLARGE=np.dot(C_mat_inv_LARGE,coord_reg_HYBRID.T)
         cif_writeALTERNATE(myid,'correctexpandedhybrid',natoms,atom_type,HYBRIDLARGE,aLEN,bLEN,cLEN,alphaD,betaD,gammaD,x_expand,y_expand,z_expand)

  ###WRITE VANILLA XYZ FILE (regular coordinates)
         returnindex3=xyz_write(myid,name,natoms,ATOM_TYPE_REAL_HYBRID,coord_reg_HYBRID)

############################################################################################################
if __name__=="__main__":
     #time.sleep(0)
     #global numprocs
     global myid
     #global node
     #numprocs = COMM.Get_size()    # Number of processors
     #myid = COMM.Get_rank()        # Id of this processor
     #node = MPI.Get_processor_name()
     #print ("I am proc %d of %d on node %s" %(myid, numprocs, node))
     #time.sleep(0)
     myid=0

     oneprint(myid,'PERIODIC BUILDER is running as an independent program.')
     oneprint(myid,'The following PERIODIC BUILDER input arugments are '+str(sys.argv[1::]))
     
     try: 
          main(myid,sys.argv[1])
     except:
          main(myid,None)
else:
     print('PERIODIC BUILDER is running as an imported module.')