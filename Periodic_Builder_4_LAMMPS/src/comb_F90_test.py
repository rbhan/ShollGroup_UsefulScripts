#This is the example script I am using to test Fortran90 modules
#Currently Loaded Modulefiles:
  #1) python/2.7(default)   2) intel/13.2.146        3) openmpi/1.5.4         4) hwloc/1.2(default)    5) mkl/11.0
#RUN COMMAND::mpirun -np 4 python comb_F90_test.py
#OUTPUT::I am proc 0 of 4 on node joe-6.pace.gatech.edu
#I am proc 1 of 4 on node joe-6.pace.gatech.edu
#I am proc 3 of 4 on node joe-6.pace.gatech.edu
#I am proc 2 of 4 on node joe-6.pace.gatech.edu
#[88889, 10, 10, 10]

from mpi4py import MPI
fcomm = MPI.COMM_WORLD.py2f() #get the Fortran communication "piece": whatever this means

import HELLLLO #get the Fortran library

numprocs = MPI.COMM_WORLD.Get_size()    # Number of processors
myid = MPI.COMM_WORLD.Get_rank()        # Id of this processor
node = MPI.Get_processor_name()
print ("I am proc %d of %d on node %s" %(myid, numprocs, node))

#Set x to 9 on all processors
x=9
if myid == 0:
  x=88888 #set x to something else only on one processor

y=None
y=HELLLLO.sayhello(fcomm,x) #run the fortran on individual processors

z=MPI.COMM_WORLD.gather(y) #gather all the "different answers"
if myid==0:
    print(z) 


#What I really want to to have an array, split it up in python, run the segements in fortran, "gather" the fortran outputs, and continue in python
#OORRRRRR I could just dump it all (the array) into the Fortran code, have it split it there, return the output already recombined to python

#I want to do whatever is easiest, efficient, and in a way that reduces the likelihood of a segmentation fault (i.e. writing to some variable/array that I shouldn't be) 