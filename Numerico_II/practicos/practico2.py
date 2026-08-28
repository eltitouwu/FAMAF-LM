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
        L[i+1:,i]=U[i+1:,i]/U[i,i]
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

print(egauss_window3(B,c),"egaus window")

#ej 10)

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

def dlup(A):
    (N,M)=A.shape
    assert(N==M)
    U=np.ndarray.copy(A)
    L=np.zeros((N,N))
    p=[i for i in range(N)]
    for i in range(N):
        for j in range(i+1,N):
            if(abs(U[p[j],i])>abs(U[p[i],i])): p[i],p[j] = p[j],p[i]
        L[p[i],i]=1
        for j in range(i+1,N):
            L[p[j],i]=U[p[j],i]/U[p[i],i]
            U[p[j],i]=0
            U[p[j],i+1:]-=U[p[i],i+1:]*L[p[j],i]
    return (L,U,p)

(R,P)=egaussp(A,b)
print(R,"egausp R de A")
print(P,"egausp P de A")

(R,P)=egaussp(B,c)
print(R,"egausp R de B")
print(P,"egausp P de B")

(L,U,P)=dlup(A)
print(L,"dlup L de A")
print(U,"dlup U de A")
print(P,"dlup P de A")

(L,U,P)=dlup(B)
print(L,"dlup L de B")
print(U,"dlup U de B")
print(P,"dlup P de B")

#ej 11)

def invalid(N):
    return np.array([np.nan]*n)

def sol_trinffilp(A,b,p):
    (N,M)=A.shape
    if(N!=b.shape[0]): return invalid(M)
    x=np.array([0]*M,dtype=np.float64)
    for i in range(M):
        if(np.abs(A[p[i]][i])<=EPS):
            if(np.abs(np.inner(A[p[i],:i],x[:i])-b[p[i]])>EPS):
                return invalid(M)
            continue
        x[i]=(b[p[i]]-np.inner(A[p[i],:i+1],x[:i+1]))/A[p[i]][i]
    return x

def sol_trinfcolp(A,_b,p):
    b=np.ndarray.copy(_b)
    (N,M)=A.shape
    if(N!=b.shape[0]): return invalid(M)
    x=np.array([0]*M,dtype=np.float64)
    for i in range(M):
        if(np.abs(A[p[i]][i])<=EPS):
            if(np.abs(b[p[i]])>EPS):
                return invalid(M)
            continue
        x[i]=b[p[i]]/A[p[i],i]
        for j in range(i+1,N)
            b[p[j]]-=A[p[j],i]*x[i]
    return x


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

def sol_trsupcolp(A,_b,p):
    b=np.ndarray.copy(_b)
    (N,M)=A.shape
    if(N!=b.shape[0]): return invalid(M)
    x=np.array([0]*M,dtype=np.float64)
    for i in range(M-1,-1,-1):
        if(np.abs(A[p[i],i])<=EPS):
            if(np.abs(b[p[i]])>EPS):
                return invalid(M)
            continue
        x[i]=b[p[i]]/A[p[i],i]
        for j in range(i):
            b[p[j]]-=A[:i,i],x[i])
    return x




def sol_egauss(_A,_b):
    (R,P)=egaussp(_A,_b)
    A=R[:,:R.shape[1]-1]
    b=R[:,R.shape[1]]
