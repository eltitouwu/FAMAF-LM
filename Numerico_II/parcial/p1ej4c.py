import numpy as np

EPS=1e-9 #epsilon de tolerancia en error absoluto




def dlup(A):
    (N,M)=A.shape
    assert(N==M)
    U=np.ndarray.copy(A)
    L=np.zeros((N,N))
    p=[i for i in range(N)]
    for i in range(N):
        for j in range(i+1,N):
            if(abs(U[p[j],i])>abs(U[p[i],i])): p[i],p[j] = p[j],p[i]
        L[p[i],p[i]]=1
        for j in range(i+1,N):
            L[p[j],p[i]]=U[p[j],i]/U[p[i],i]
            U[p[j],i]=0
            U[p[j],i+1:]-=U[p[i],i+1:]*L[p[j],p[i]]
    return (L,U,p)


def invalid(n):
    return np.array([np.nan]*n)


def sol_trsupfilp(A,b,p):
    (N,M)=A.shape
    if(N!=b.shape[0]): return invalid(M)
    x=np.array([0]*M,dtype=np.float128)
    for i in range(M-1,-1,-1):
        if(np.abs(A[p[i]][i])<=EPS):
            if(np.abs(np.inner(A[p[i],i+1:],x[i+1:])-b[p[i]])>EPS):
                return invalid(M)
            continue
        x[i]=(b[p[i]]-np.inner(A[p[i],i+1:],x[i+1:]))/A[p[i],i]
    return x

def sol_trinffilpp(A,b,p):
    (N,M)=A.shape
    if(N!=b.shape[0]): return invalid(M)
    x=np.array([0]*M,dtype=np.float128)
    for i in range(M):
        aux=0.
        for j in range(i):
            aux+=A[p[i],p[j]]*x[p[j]]
        
        if(np.abs(A[p[i]][p[i]])<=EPS):
            if(np.abs(aux-b[p[i]])>EPS):
                return invalid(M)
            continue

        x[p[i]]=(b[p[i]]-aux)/A[p[i]][p[i]]
    return x




def inv_Ak_kLU(A,k):
    (L,U,p)=dlup(A)
    X=np.eye(A.shape[0])
    for i in range(A.shape[0]):
        for t in range(k):
            X[:,i]=sol_trsupfilp(U,sol_trinffilpp(L,X[:,i],p),p)
    return X

A=np.eye(5,dtype=np.float128)*2
for i in range(4):
    A[i+1,i]=-1
    A[i,i+1]=-1

X=inv_Ak_kLU(A,3)


assert(np.allclose(A@A@A@X,np.eye(5),atol=EPS))

