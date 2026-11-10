import numpy as np

EPS=1e-9 #epsilon de tolerancia en error absoluto

def invalid(n):
    return np.array([np.nan]*n)

def sol_trinffil(A,b):
    (n,m)=A.shape
    if(n!=b.shape[0]): return invalid(m)
    x=np.array([0]*m,dtype=np.float128)
    for i in range(m):
        if(np.abs(A[i][i])<=EPS):
            if(np.abs(np.inner(A[i,:i],x[:i])-b[i])>EPS):
                return invalid(m)
            continue
        x[i]=(b[i]-np.inner(A[i,:i+1],x[:i+1]))/A[i][i]
    return x

def sol_trsupfil(A,b):
    (n,m)=A.shape
    if(n!=b.shape[0]): return invalid(m)
    x=np.array([0]*m,dtype=np.float128)
    p=min(n,m)
    for i in range(p-1,-1,-1):
        if(np.abs(A[i][i])<=EPS):
            if(np.abs(np.inner(A[i,i+1:],x[i+1:])-b[i])>EPS):
                return invalid(m)
            continue
        x[i]=(b[i]-np.inner(A[i,i+1:],x[i+1:]))/A[i][i]
    return x

def cholesky(_X):
    X=np.ndarray.copy(_X)
    n=X.shape[0]
    assert(n==X.shape[1])
    G=np.zeros((n,n),dtype=np.float128)
    for i in range(n):
        assert(X[i][i]>EPS)
        G[i][i]=np.sqrt(X[i][i])
        G[i,i+1:]=X[i,i+1:]/G[i][i]
        X[i+1:,i+1:]-=np.outer(G[i,i+1:],G[i,i+1:])
    return G



def solve_ej1(C,B,b):
    G1=cholesky(C)
    
    K=np.zeros(B.shape,dtype=np.float128)
    for i in range(B.shape[1]):
        K[:,i]=sol_trinffil(G1.T,B[:,i])
    G2=cholesky(K.T@K)
    Izq=np.block([[G1.T,np.zeros(B.shape)],[K.T,G2.T]])
    Der=np.block([[G1,K],[np.zeros((B.T).shape),-G2]])

    return sol_trsupfil(Der,sol_trinffil(Izq,b))

C=np.array([[1,2,1],[2,5,2],[1,2,5]],dtype=np.float128)
B=np.array([[1,0],[4,3],[1,8]],dtype=np.float128)
b=np.array([5,16,17,6,11],dtype=np.float128)
x=solve_ej1(C,B,b)
A=np.block([[C,B],[B.T,np.zeros((2,2))]])

# print(A@x,b)

assert(np.allclose(A@x,b,rtol=EPS))

