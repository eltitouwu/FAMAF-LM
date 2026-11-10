import numpy as np
from sol_egauss import sol_trsupfilp
from sol_tr import invalid
EPS=np.float128(10)**-20

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

def sol_dlu(X,b):
    (L,U,P)=dlup(X)
    # assert(np.allclose(L@U,A,atol=EPS,rtol=EPS))
    return sol_trsupfilp(U,sol_trinffilpp(L,b,P),P)