import numpy as np
EPS=np.float128(10)**-20
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

def sol_trinfcol(A,_b):
    b=np.ndarray.copy(_b)
    (n,m)=A.shape
    if(n!=b.shape[0]): return invalid(m)
    x=np.array([0]*m,dtype=np.float128)
    for i in range(m):
        if(np.abs(A[i][i])<=EPS):
            if(np.abs(b[i])>EPS):
                return invalid(m)
            continue
        x[i]=b[i]/A[i][i]
        b[i+1:]-=np.inner(A[i+1:,i],x[i])
    return x

#b)

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

def sol_trsupcol(A,_b):
    b=np.ndarray.copy(_b)
    (n,m)=A.shape
    if(n!=b.shape[0]): return invalid(m)
    x=np.array([0]*m,dtype=np.float128)
    for i in range(m-1,-1,-1):
        if(np.abs(A[i][i])<=EPS):
            if(np.abs(b[i])>EPS):
                return invalid(m)
            continue
        x[i]=b[i]/A[i][i]
        b[:i]-=np.inner(A[:i,i],x[i])
    return x

def cholesky(_X):
    X=np.ndarray.copy(_X)
    n=X.shape[0]
    assert(n==X.shape[1])
    G=np.zeros((n,n),dtype=X.dtype)
    for i in range(n):
        assert(X[i][i]>EPS)
        G[i][i]=np.sqrt(X[i][i])
        G[i,i+1:]=X[i,i+1:]/G[i][i]
        X[i+1:,i+1:]-=np.outer(G[i,i+1:],G[i,i+1:])
    return G

def sol_defpos(X,b):
    G=cholesky(X)
    return sol_trsupfil(G,sol_trinffil(G.T,b))