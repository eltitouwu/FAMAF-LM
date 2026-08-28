import numpy as np

#ej 6
m=9
n=6
A=np.random.rand(m,n)

#a)

#i)
print(A[:,n-1],A[m-1,:])

#ii)

print(A[:m-2,n-1],A[m-2,:n-2])

#iii)

print(A[n-2:m-2,m-n-1])

#iv)

print(A[n:m,2:n])


#b) 
def generador_de_matrices_en_bloques(X):
    (N,M)=X.shape()
    R=[]
    for i in range(min(N,M)):
        R.append(X[:i+1,:i+1])    
    return R

#ej 7


A=np.array([[2,-1,0],[-1,2,-1],[0,-1,2]])
B=np.array([[4,-1,0],[-1,4,-1],[0,-1,4]])

#a)

C=np.block([[A,-np.eye(3),np.zeros((3,3))],[-np.eye(3),B,-np.eye(3)],[np.zeros((3,3)),-np.eye(3),A]])


#b)
X1=np.block([[np.eye(3),np.zeros((3,3)),np.zeros((3,3))],[np.zeros((3,3)),np.eye(3),np.zeros((3,3))]]);
X1=X1.T


X2=np.block([np.zeros((3,3)),np.eye(3),np.zeros((3,3))]);
X2=X2.T


X3=np.block([[np.zeros((3,3)),np.eye(3),np.zeros((3,3))],[np.zeros((3,3)),np.zeros((3,3)),np.eye(3)]])
X3=X3.T

#c)
assert(np.array_equal(X1.T@C@X1,C[:6,:6]))
assert(np.array_equal(X2.T@C@X2,C[3:6,3:6]))
assert(np.array_equal(X3.T@C@X3,C[3:,3:]))


