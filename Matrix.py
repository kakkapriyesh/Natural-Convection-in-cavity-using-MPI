from mpi4py import MPI
import numpy as np
i=4;
j=4;
k=2;
global data3;
global data4;
comm=MPI.COMM_WORLD 
rank=comm.Get_rank()
n=comm.Get_size()
n=2
if (rank==0):
 A=np.random.rand(i,k)
 F=[]
 print ("Matrix F is \n",F)
 #B=np.ones(shape=(k,j))
 B=np.random.rand(k,j)
 #print "Matrix B is \n",B
 r=len(A)
 c=len(B[0])
 print ("the row of C is ",r,"the columns of C is ",c)
 P=np.dot(A,B), "Big matrix"
 print (P)
 OR=float(r)/c;
 #print ("Aspect ratio for matrix is ",OR;)
 diff=1;
 diff2=23;
 C="hi"
 #mapping nprocs ro C
 for I in range (1,n):
   
  if n%I==0:
   

    AR=float(I*I)/n
   
    diff = abs(float(AR-OR))
    if (diff<diff2):
         diff2=diff
    print (diff ,"AR ", AR,"iproc ",I,"jproc ",(n/I))
    rn= I
    cn= int (n/I) 
 K=int (r/rn)
 J=int (c/cn)

 print ("size of row of A matrix is ",K,"size of coloumn of B matrix is ",J)
 D=np.linspace(0,j,(cn+1),dtype=int)
 #print D
 E=np.linspace(0,i,(rn+1),dtype=int)
 print (n,"fgh",rn)

# if (rank==0):
 for m in range(0,rn):
   for n in range(0,cn):
    data1=np.array(A[E[m]:E[m+1],:k])
    comm.send([data1],dest=(m*cn+n+1),tag=11)
    data2=np.array(B[:k,D[n]:D[n+1]])
    comm.send([data2],dest=(m*cn+n+1),tag=12)
    #C=0
    #print data1 ," A1"
    #print data2, " B1"
elif (rank!=0):
     #data1=np.empty([1],dtype='float32');
     data3=comm.recv(source=0,tag=11);
     #data2=np.empty([1],dtype='float32');
     data4=comm.recv(source=0,tag=12);
     C=np.dot(data3,data4)
     
     print((C.shape),"see")
     print (C," for rank ", rank)
     
matrix=comm.gather(C,root=0)
np.delete(matrix,[0])
print (matrix,"sthis is the matrix")
    if rank == 0:
		C = np.empty([i, j], dtype = 'float32')
		k = 0
		for i in range(int(rn):
			for j in range(int(cn)):
				C[E[i]:E[i+1], D[j]:D[j+1]] = matrix[k]
				k += 1
		print("Shape of Matrix C (Matrix A * Matrix B) is : {}".format(C.shape))
		print ""
		print("Matrix C is:")
		print(C)
    
    
    
    
    
    
    
    
    
    
    
    
    
    # comm.send(matrix,dest=0,tag=13)
#L=np.array(L)
    # print(L)
#K=np.reshape(L,(i,j))

#data = np.dot(data3, data4)
	#print(data.shape, rank)

#matrix = comm.gather(data, root = 0)
	# To print final matrix in rank 0
'''if (rank == 0):
     data=comm.recv(source=0,tag=13)
     #matrix = comm.gather(data, root = 0)
     print (data)'''
		
