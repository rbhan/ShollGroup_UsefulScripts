from collections import deque
import numpy as np
from collections import defaultdict
import time

def linkerconnection_format(METAL_DICT_HYBRID):
#Get the NN list for the Zn atoms
  ZnNNdict=defaultdict(list)
  ZnKEYS=list()
  for clust in range(len(METAL_DICT_HYBRID)):
    Znbond=METAL_DICT_HYBRID[clust+1]
    #print(Znbond)
    ZnNNdict[Znbond[0]].append(Znbond[1])
    ZnNNdict[Znbond[1]].append(Znbond[0])
    if Znbond[0] not in ZnKEYS:
      ZnKEYS.append(Znbond[0])
    if Znbond[1] not in ZnKEYS:
      ZnKEYS.append(Znbond[1])
  ZnKEYS=np.sort(ZnKEYS)  

  #print(len(ZnNNdict))
  #print(len(ZnKEYS))

  molecule=[]
  for Zn in range(len(ZnNNdict)):
    dicttemp={}
    dicttemp["connections"]=list(ZnNNdict[ZnKEYS[Zn]])
    dicttemp["key"]=ZnKEYS[Zn]
    molecule.append(dicttemp)

#pentagon = [{"key":"A", "connections":["B","E"]},
#   {"key":"B", "connections": ["A", "C"]},
#   {"key":"C", "connections": ["B", "D"]},
#   {"key":"D", "connections": ["C", "E"]},
#   {"key":"E", "connections": ["D", "A"]}]

  #print(molecule)
  return(molecule,ZnNNdict)

def postprocessrings(rings,METAL_DICT_HYBRID,ZnNNdict,COM_cluster_TYPES):
  ringdict=defaultdict(list)
  ringcounter=0
  for ring in rings:
    #print('ring',ring)
    ZNbond1=[]
    ZNbond2=[]
    for element in ring:
      #print('element',element)
      for Zn in ZnNNdict[element]: #access that Zn atoms NNs
        #print('Zn',Zn)
        if Zn in ring:
          #print('attachedZn',Zn)
          ZNbond1.append(element)
          ZNbond2.append(Zn)
    ringdict[ring].append(ZNbond1)
    ringdict[ring].append(ZNbond2)
    #print(ringdict)
    
  #print(ringdict)

  #print(len(ringdict)) #this dictionary now has the list of bonds, duplicates of which can be removed

  RINGTYPES=defaultdict(list)
  for ring in rings:
    listoflists=ringdict[ring]
    #print(listoflists)
    bond1=listoflists[0]
    bond2=listoflists[1]
    stackedbonds=np.vstack((bond1,bond2))
    stackedbonds=np.transpose(stackedbonds)
    stackedbonds=np.sort(stackedbonds)
    #print(stackedbonds)
    a=stackedbonds
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)
    unique_a = a[idx]
    #print(unique_a)

    
    #ORDER THE BONDS IN THE CORRECT CONNECTIVITY
    unique_aB=[]
    rowcounter=0
    for row1 in unique_a:
      if rowcounter==0:
        unique_aB.append(row1)
      if rowcounter==5: #----------------------------->possible input parameter
         break

      for row2 in unique_a:
        #print('ROW',rowcounter,row1,row2)

        CURRENTROW=unique_aB[rowcounter]
        #print(CURRENTROW,row2)

        if all(x in CURRENTROW for x in row2):
          BLAH='NONE'
        else:
            #print('PASSED',CURRENTROW,row2)
            if CURRENTROW[1]==row2[0]: #look at the end of the array
              unique_aB.append(row2)
            elif CURRENTROW[1]==row2[1]: #reverse the array in appending
              unique_aB.append(row2[::-1])

      rowcounter+=1
    #print('ARRANGED',unique_aB)

    #THIS LABELS THE "BONDS"    
    for row in range(len(unique_aB)):
      actualbondrow=list(unique_aB[row])
      for clust in range(len(METAL_DICT_HYBRID)):
        Znbond=list(METAL_DICT_HYBRID[clust+1])
        if all(x in Znbond for x in actualbondrow):
          CONNECTIONTYPE=(COM_cluster_TYPES[clust+1])
          CONNECTIONTYPE=int(CONNECTIONTYPE[0])
          #print(Znbond,actualbondrow,CONNECTIONTYPE)
          RINGTYPES[ring].append(CONNECTIONTYPE)

    #With unique a I can determine the number of a certain type of linker

  print(RINGTYPES) #THIS DOES NOT MAINTAIN ORDER
  return(RINGTYPES)

def histo(rings,RINGTYPES):
  HISTOGRAMCOUNTS=[]
  for ring in rings:
    linkers=RINGTYPES[ring]
    sum=np.sum(linkers)
    HISTOGRAMCOUNTS.append(sum)
  HISTOGRAMBIN=np.bincount(HISTOGRAMCOUNTS)
  return(HISTOGRAMBIN,HISTOGRAMCOUNTS)

         
def getring(Molecule, rootnode):
    # connections is just the table of connections in the molecule
    connections = { mol['key']: mol['connections'] for mol in Molecule}
    # ring members must be a member of the ring
    assert(len(connections[rootnode])>1)
    # paths is a dictionary of paths starting form the start node. It is initialized
    # to zero and will be used to updated as each member of the queue is analyzed
    paths = {mol['key']: set() for mol in Molecule}
    
    # atomqueue will keep track of atoms that have been visited and need to be processed
    atomqueue = deque()
    
    # visited is the atom keys that have been seen already
    visited = set()
    
    #when a ring is found, return the ring
    ring = None
    
    # initialize the root node
    atomqueue.append(rootnode)
    paths[rootnode].add(rootnode)
    visited.add(rootnode)

    while not ring:
        #print "atomqueue: {}".format(atomqueue)
        #print "visited: {}".format(visited)
        #print "paths: {}".format(paths)
        
        try:
            atomcode = atomqueue.popleft()
            #print "atomcode: {}".format(atomcode)

            connectedatoms = [a for a in connections[atomcode] if a not in visited]
            #print "atoms connected: {}".format(connectedatoms)
            for atom in connectedatoms:
                #print "atom: {}".format(atom)
                #check for rings and return if a ring is found
                currentpath = paths[atomcode]
                nextpath  = paths[atom]
                intersect = currentpath.intersection(nextpath)
                #print "current path: {}".format(currentpath)
                #print "next path: {}".format(currentpath)
                if len(nextpath) > 3 and len(intersect) == 1:
                    return currentpath.union(nextpath)  ##############this is the final statement
    
                # update path only if the target path is empty 
                # to avoid incorporating longer paths
                if len(nextpath) == 0:
                    paths[atom] = currentpath.copy()
                    paths[atom].add(atom)
            
                #add atoms to queue
                atomqueue.append(atom)
                
            #update visited
            visited.add(atomcode)
        except:
            raise Error("Problem Here. You should never get here!")

def getrings(Molecule):
    """
    find atoms that are connected by only a single bond,
    eliminate them and the connections to them by other atoms until
    there are only atoms that have at least two connections
    then find the rings on them
    """
    
    # copy the Molecule for local editing
    tempmol = list(Molecule)
    #print tempmol
    # find the atoms that are initially less than 2 bond connections
    singles_index = [i for i,a in enumerate(tempmol) if len(a['connections']) < 2]
    
    
    # eliminate all atoms that are not part of a ring
    while len(singles_index) > 0:
        for idx in singles_index:
            key = tempmol[idx]['key']
            connections = tempmol[idx]['connections']
            
            # remove connections 
            for con in connections:
                #print "index: {}\nkey:{}\nconnection:{}".format(idx,key,con)
                index = [i for i,x in enumerate(tempmol) if x['key'] == con][0]
                tempmol[index]['connections'].remove(key)
 
            # remove atom
            tempmol.pop(idx) 
            #print tempmol
            
            # reset singles index
            singles_index = [i for i,a in enumerate(tempmol) if len(a['connections']) < 2]
            #print len(singles_index)
        
    # calculate ring on the trimmed molceule
    rings = set()
    for atom in tempmol:
        if len(atom['connections']) > 1:
            ring = frozenset(getring(Molecule, atom["key"]))
            #print atom, ring
            if ring not in rings:
                rings.add(ring)
    return(rings)

