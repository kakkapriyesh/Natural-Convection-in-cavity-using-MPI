from mpi4py import MPI
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
import numpy as np
n=comm.Get_size();


 
if(rank==0):
     l=int(np.random.random(1)*10)
     a=np.arange(l)
     print (l,": this is the data size")
     print (a,": this is the array")     
     data2=np.array(l,dtype='i')
     data1=np.array(a,dtype='i')
     comm.Send([data2,MPI.INT],dest=1,tag=11)
     comm.Send([data1,MPI.INT],dest=1,tag=11)
elif(rank==1):
    data2=np.empty([1],dtype='i')
    comm.Recv([data2,MPI.INT],source=0,tag=11)
    data1=np.empty([data2[0]],dtype='i')
    comm.Recv([data1,MPI.INT],source=0,tag=11)
    print (data1,"recived at ",rank,"from ",0 ) 
 
