
import numpy as np
from numpy import matlib as mb
import sys
def gcd(a, b):
    # Return the GCD of a and b using Euclid's Algorithm
    while a != 0:
        a, b = b % a, a
    return b


def findModInverse(a, m):
    # Returns the modular inverse of a % m, which is
    # the number x such that a*x % m = 1

    if gcd(a, m) != 1:
        return None # no mod inverse if a & m aren't relatively prime

    # Calculate using the Extended Euclidean Algorithm:
    u1, u2, u3 = 1, 0, a
    v1, v2, v3 = 0, 1, m
    while v3 != 0:
        q = u3 // v3 # // is the integer division operator
        v1, v2, v3, u1, u2, u3 = (u1 - q * v1), (u2 - q * v2), (u3 - q * v3), v1, v2, v3
    return u1 % m

def codingB(B:np.ndarray,p:int,n:int,m:int,x:int):
    (row,col)=B.shape
    a=row//p
    b=col//n
    B_i=mb.zeros((a,b))
    for j in range(p):
        for k in range(n):
            # print(A[j*a:j*a+a,k*b:k*b+b])
            B_i=B_i+np.multiply(B[j*a:j*a+a,k*b:k*b+b],(x**(p-1-j+k*p*m))%10007)   
            np.mod(B_i,10007)
    print("B"+str(x)+"=",B_i)
    return B_i

def codingA(A:np.ndarray,p:int,m:int,x:int):
    (row,col)=A.shape
    a=row//p
    b=col//m
    A_i=mb.zeros((a,b))
    for j in range(p):
        for k in range(m):
            # print(A[j*a:j*a+a,k*b:k*b+b])
            A_i=A_i+np.multiply(A[j*a:j*a+a,k*b:k*b+b],(x**(j+k*p))%10007)
            np.mod(A_i,10007)
    print("A"+str(x)+"=",A_i)
    return A_i


def getCiRowCol(A,B,p,m,n):
    a=A.shape[1]//m
    b=B.shape[1]//n
    return (a,b)
    
def getpoly(x_i,x_j):
    inv=findModInverse((x_i-x_j+10007)%10007,10007)
    if inv is not None:
        return np.poly1d(np.array([(inv+10007)%10007,(-inv*x_j+10007)%10007]))
    else:
        return None

def decoding(_C,x,p,m,n):
    (a,b)=_C[0].shape
    coeM=[]
    for i in range(p*m*n+p-1) :
        coeM.append(np.empty([a,b]))
    print("row:",a,"col:",b)
    for c_i in range(a):
        for c_j in range(b):
            f=np.poly1d(np.array([0])) 
            for i in range(p*m*n+p-1):
                subf=np.poly1d(np.array([1])) 
                for j in range(p*m*n+p-1):
                    
                    if j!=i:
                        subsubf=getpoly(x[i],x[j])
                        if subsubf is not None:
                            subf=np.convolve(subf,subsubf)
                            subf=np.mod(subf,10007)
                        else:
                            print("===============error================")
                            sys.exit()
                f=f+np.mod(subf*_C[i][c_i,c_j],10007)
                f=np.mod(f,10007)
            # print(f)
            for k in range(p*m*n+p-1):
                coeM[k][c_i,c_j]=f[k]
    for i in range(m):
        for j in range(n):
            print("C("+str(i)+","+str(j)+")=")
            print(coeM[p-1+(m-i-1)*p+(n-j-1)*p*m])
    # print(coeM)



print("===============input===============")
A = np.array([[1, 2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]])
B = np.array([[1,2,3,4,5,6],[7,8,9,10,11,12],[13,14,15,16,17,18],[19,20,21,22,23,24]])
print("A=")
print(A)
print("B=")
print(B)

print("\n===============(A^T)*B:===============")
print (np.dot(A.T,B))


p=2
m=2
n=3
_A:list[np.ndarray]=[]
_B:list[np.ndarray]=[]
_C:list[np.ndarray]=[]
x=[]
k=p*m*n+p-1
print("\n===============coding===============")
for i in range(1,k+1):
    x.append(i)
    _A.append(codingA(A,p,m,i)) 
    _B.append(codingB(B,p,n,m,i)) 
    _C.append(np.dot(_A[i-1].T,_B[i-1]))
    print("C"+str(i)+"=",_C[i-1])

print("\n=================decoding==============")
decoding(_C,x,p,m,n)


