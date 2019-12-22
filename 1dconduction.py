from mpi4py import MPI
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
import numpy as np
nproc=comm.Get_size();
#l=500.0;
n=6;
x=1.0;
print 
U=np.zeros(n)
Ut=np.zeros(n)
print rank;
print U;
t=1;
time=0;
U[0]=0
U[n-1]=1
Ut[n-1]=1
d=(n/(nproc))+2;
U=np.zeros(d)
print d
while time<10000:

  

  
   
  if(rank==0):
     
     for i in range (1,d-1,1):
       U=np.zeros(d-1)
       U[0]=1;
       print (U,"fgfw")
       
       print "sfds"
       Ut[i]=U[i]+t*0.5*(U[i+1]-2*U[i]+U[i-1])/(x*x)
       print (U,"fgf")
     U=Ut[1:d]
     data1=np.empty([1],dtype='float32')
     comm.Recv([data1,MPI.INT],source=1,tag=11)
     U[d-1]=data1
     data=np.array(U[d-2],dtype='float32')
     comm.Send([data,MPI.INT],dest=1,tag=11)
     print "sfds1"
  elif(rank==1):
     for i in range (0,d-1,1):
      U=np.zeros(d-1)
      U[d-1]=1

      print (U,"sfds10")     
      Ut[i]=U[i]+t*0.5*(U[i+1]-2*U[i]+U[i-1])/(x*x)
     U=Ut[d:2*d]     
     data=np.empty([1],dtype='i')
     comm.Recv([data,MPI.INT],source=0,tag=11)
     U[d-1]=data
     print "sfds2"
     data1=np.array(U[d],dtype='i')
     comm.Send([data1,MPI.INT],dest=1,tag=11)
     print (data,"recived at ",rank,"from ",0 ) 
 



   
     U=Ut[:]
  time=time+t
print (U)




