import numpy as np

EPS=1e-9

#ej 5)

def egauss(_A,_b):
    (N,M)=_A.shape
    assert(N==_b.shape[0])
    A=np.block([np.ndarray.copy(_A),_b.reshape(-1,1)])
    for i in range(N):
        assert(abs(A[i,i])>EPS)
        A[i,i+1:]/=A[i,i]
        A[i,i]=1
        A[i+1:,i+1:]-=np.outer(A[i+1:,i],A[i,i+1:])
        A[i+1:,i]=np.zeros((1,N-i-1),dtype=np.float64)
    return A


A=np.array([[1,5,7],[2,8,3],[1,1,1]],dtype=np.float64)
b=np.array([1,2,3],dtype=np.float64)
B=np.array([[1,5,0,7],[2,8,3,0],[0,1,1,1],[9,0,6,7]],dtype=np.float64)
c=np.array([1,2,3,4],dtype=np.float64)

print(A,"A")
print(b,"b")
print(egauss(A,b),"egaus")
print(B,"B")
print(c,"c")
print(egauss(B,c),"egauss")


#ej 6)

def dlu(A):
    (N,M)=A.shape
    assert(N==M)
    U=np.ndarray.copy(A)
    L=np.eye(N)
    for i in range(N):
        assert(abs(U[i,i])>EPS)
        L[i+1:,i]=U[i+1:,i]/U[i][i]
        U[i+1:,i]=np.zeros((1,N-i-1),dtype=np.float64)
        U[i+1:,i+1:]-=np.outer(L[i+1:,i],U[i,i+1:])
    return (L,U)

(L,U)=dlu(A)
print(L,"L de A")
print(U,"U de A")
assert(np.allclose(L@U,A,atol=EPS,rtol=EPS))
(L,U)=dlu(B)
print(L,"L de B")
print(U,"U de B")
assert(np.allclose(L@U,B,atol=EPS,rtol=EPS))
#ej 7)

def egauss_window3(_A,_b):
    (N,M)=_A.shape
    assert(N==_b.shape[0] and N==M)
    A=np.block([np.ndarray.copy(_A),_b.reshape(-1,1)])
    for i in range(N):
        assert(abs(A[i,i])>EPS)
        id=[i+1]
        if(i+1<N-1): id.append(N-1)
        if(i+1<N): id.append(N)
        for j in id:
            A[i,j]/=A[i,i]
        A[i,i]=1
        for j in id:
            if(j>=N): break
            for k in id:
                A[j,k]-=A[i,k]*A[j,i]
            A[j,i]=0
    A[N-1,N]/=A[N-1,N-1]
    A[N-1,N-1]=1
    return A

print(egauss_window3(A,b),"egauss window")

A=np.array([[1,5,7],[2,8,3],[1,1,1]],dtype=np.float64)
b=np.array([1,2,3],dtype=np.float64)

print(egauss_window3(B,c),"egaus window")