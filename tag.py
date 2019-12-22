from mpi4py import MPI
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
import numpy as np
count=np.array([0],dtype='i')
n=comm.Get_size();
if (rank==0):
  data1=np.array(['35','JOHN','67.8'],dtype='str')
  #data2=np.array([2],dtype='i')
  #a=((rank+1)%n);
  comm.Send([data1,MPI.CHAR],dest=1,tag=11)
  #comm.Send([data2,MPI.INT],dest=a,tag=10)
  #print (data1,"and",data2,"sent from ",rank,"to ",rank+1 )
elif(rank==1):

  data1=np.empty([1],dtype='str')
  #data2=np.empty([1],dtype='i')
  comm.Recv([data1,MPI.CHAR],source=0,tag=11)
  #comm.Recv([data2,MPI.CHAR],source=(rank-1)%n,tag=10)
  print (data1,"recived at ",rank,"from ",rank-1 )
  
 
