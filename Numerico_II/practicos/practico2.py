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
        L[p[i],p[i]]=1
        for j in range(i+1,N):
            L[p[j],p[i]]=U[p[j],i]/U[p[i],i]
            U[p[j],i]=0
            U[p[j],i+1:]-=U[p[i],i+1:]*L[p[j],p[i]]
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
    return np.array([np.nan]*N)

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
        for j in range(i+1,N):
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
            b[p[j]]-=A[p[j],i]*x[i]
    return x




def sol_egauss(_A,_b):
    (R,P)=egaussp(_A,_b)
    A=R[:,:R.shape[1]-1]
    b=R[:,R.shape[1]-1]
    return sol_trsupfilp(A,b,P);

A=np.array([[2,10,8,8,6],[1,4,-2,4,-1],[0,2,3,2,1],[3,8,3,10,9],[1,4,1,2,1]],dtype=np.float64);
b1=np.array([52,14,12,51,15],dtype=np.float64)
b2=np.array([50,4,12,48,12],dtype=np.float64)
assert(np.allclose(A@sol_egauss(A,b1),b1,atol=EPS,rtol=EPS))
assert(np.allclose(A@sol_egauss(A,b2),b2,atol=EPS,rtol=EPS))

#12)

def sol_trinffilpp(A,b,p):
    (N,M)=A.shape
    if(N!=b.shape[0]): return invalid(M)
    x=np.array([0]*M,dtype=np.float64)
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

def inv_lu(X):
    (L,U,P)=dlup(X)
    # assert(np.allclose(L@U,A,atol=EPS,rtol=EPS))
    Y=np.zeros(X.shape)
    ei=np.zeros(X.shape[0])
    for i in range(X.shape[0]):
        ei[i]=1
        Y[:,i]=sol_trsupfilp(U,sol_trinffilpp(L,ei,P),P)
        ei[i]=0
    return Y


assert(np.allclose(A@inv_lu(A),np.eye(A.shape[0]),atol=EPS,rtol=EPS))


def det_lu(X):
    (N,M)=A.shape
    assert(N==M)
    U=np.ndarray.copy(A)
    L=np.zeros((N,N))
    p=[i for i in range(N)]
    det=1.
    for i in range(N):
        for j in range(i+1,N):
            if(abs(U[p[j],i])>abs(U[p[i],i])): p[i],p[j] = p[j],p[i]
        if(abs(U[p[i],i])<=EPS): return 0.
        det*=U[p[i],i]
        L[p[i],p[i]]=1
        for j in range(i+1,N):
            L[p[j],p[i]]=U[p[j],i]/U[p[i],i]
            U[p[j],i]=0
            U[p[j],i+1:]-=U[p[i],i+1:]*L[p[j],p[i]]
    return det

import time
mitime=time.time()
midet=det_lu(A)
mitime=time.time()-mitime

sutime=time.time()
sudet=np.linalg.det(A)
sutime=time.time()-sutime

assert(abs(midet-sudet)<=EPS)

print("det_lu tadó "+str(mitime)+" y np.linalg.det tardó "+str(sutime))


import matplotlib.pyplot as plt #por favor importar para graficar




def ecuacion_esfera(p1,p2,p3,p4):
    A=np.block([[p1,np.eye(1)],[p2,np.eye(1)],[p3,np.eye(1)],[p4,np.eye(1)]])
    b=-np.array([np.inner(p1,p1),np.inner(p2,p2),np.inner(p3,p3),np.inner(p4,p4)])
    v=sol_egauss(A,b)

    c=-np.array([v[0],v[1],v[2]])/2
    r=np.sqrt(c*c-v[3])

    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    N=100
    h=2*np.pi/(N)
    th=[i*h for i in range(N)]
    h/=2
    ph=[i*h for i in range(N)]
    x=[r*np.sin(b)*np.cos(a)+c[0] for a in th for b in ph]
    y=[r*np.sin(b)*np.sin(a)+c[1] for a in th for b in ph] #eje y.
    z=[r*np.cos(b)+c[2] for a in th for b in ph]
    ax.scatter3D(x,y,z) #grafica
    plt.show()    #muestro la grafica

    return v

p1=np.array([1,1,0])
p2=np.array([3,-1,3])
p3=np.array([-1,3,-4])
p4=np.array([14,84,58])

ecuacion_esfera(p1,p2,p3,p4)
