import numpy as np

EPS=1e-9

def egaussp(_A,_b):
    (N,M)=_A.shape
    assert(N==_b.shape[0])
    A=np.block([np.ndarray.copy(_A),_b.reshape(-1,1)])
    p=[i for i in range(N)]
    for i in range(N):
        for j in range(i+1,N):
            if(abs(A[p[j],i])>abs(A[p[i],i])): p[i],p[j] = p[j],p[i]
        if(abs(A[p[i],i])<=EPS): continue
        A[p[i],i+1:]/=A[p[i],i]
        A[p[i],i]=1
        for j in range(i+1,N):
            A[p[j],i+1:]-=A[p[i],i+1:]*A[p[j],i]
            A[p[j],i]=0
    return (A,p)

def invalid(N):
    return np.array([np.nan]*N)

def sol_trsupfilp(A,b,p):
    (N,M)=A.shape
    if(N!=b.shape[0]): return invalid(M)
    x=np.array([0]*M,dtype=np.float64)
    for i in range(M-1,-1,-1):
        if(np.abs(A[p[i]][i])<=EPS):
            if(np.abs(np.inner(A[p[i],i+1:],x[i+1:])-b[p[i]])>EPS):
                return invalid(M)
            continue
        x[i]=(b[p[i]]-np.inner(A[p[i],i+1:],x[i+1:]))/A[p[i],i]
    return x


def sol_egauss(_A,_b):
    (R,P)=egaussp(_A,_b)
    A=R[:,:R.shape[1]-1]
    b=R[:,R.shape[1]-1]
    return sol_trsupfilp(A,b,P);