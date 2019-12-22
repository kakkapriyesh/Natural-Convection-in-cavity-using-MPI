from mpi4py import MPI
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
import numpy as np
count=np.array([0],dtype='i')
n=comm.Get_size();
if (rank==0):
  a=((rank+1)%n);
  comm.Send([count,MPI.INT],dest=a,tag=11)
elif(rank!=0):
  comm.Recv([count,MPI.INT],source=(rank-1)%n,tag=11)
  print ("data recived from ",(rank-1)%n,"by ",rank,"and sent to ",(rank+1)%n )
  comm.Send([count,MPI.INT],dest=(rank+1)%n,tag=11)
 
